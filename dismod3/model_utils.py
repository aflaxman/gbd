import pylab as pl
import numpy as np
import pymc as mc
from pymc import gp

from dismod3.settings import *
from bayesian_models import probabilistic_utils
from bayesian_models.probabilistic_utils import trim, uninformative_prior_gp, NEARLY_ZERO, MAX_AGE

MISSING = -99

def indices_for_range(age_mesh, age_start, age_end):
    return [ ii for ii, a in enumerate(age_mesh) if a >= age_start and a <= age_end ]

def generate_prior_potentials(prior_str, age_mesh, rate, confidence_stoch):
    """
    return a list of potentials that model priors on the rate_stoch

    prior_str may have lines in the following format:
      smooth <tau> <age_start> <age_end>
      zero <age_start> <age_end>
      confidence <mean> <tau>
      increasing <age_start> <age_end>
      decreasing <age_start> <age_end>
      convex_up <age_start> <age_end>
      convex_down <age_start> <age_end>
      unimodal <age_start> <age_end>
    
    for example: 'smooth .1 \n zero 0 5 \n zero 95 100'

    age_mesh[i] indicates what age the value of rate[i] corresponds to

    confidence_stoch can be an additional stochastic variable that is used
    in the beta-binomial model, but it is not required
    """

    def derivative_sign_prior(rate, prior, deriv, sign):
        age_start = int(prior[1])
        age_end = int(prior[2])
        age_indices = indices_for_range(age_mesh, age_start, age_end)
        @mc.potential(name='deriv_sign_{%d,%d,%d,%d}^%s' % (deriv, sign, age_start, age_end, rate))
        def deriv_sign_rate(f=rate,
                            age_indices=age_indices,
                            tau=1000.,
                            deriv=deriv, sign=sign):
            df = np.diff(f[age_indices], deriv)
            return -tau * np.dot(df**2, (sign * df < 0))
        return [deriv_sign_rate]

    priors = []
    for line in prior_str.split('\n'):
        prior = line.strip().split()
        if len(prior) == 0:
            continue
        if prior[0] == 'smooth':
            tau_smooth_rate = float(prior[1])

            if len(prior) == 4:
                age_start = int(prior[2])
                age_end = int(prior[3])
            else:
                age_start = 0
                age_end = MAX_AGE
            age_indices = indices_for_range(age_mesh, age_start, age_end)
                
            @mc.potential(name='smooth_{%d,%d}^%s' % (age_start, age_end, rate))
            def smooth_rate(f=rate, age_indices=age_indices, tau=tau_smooth_rate):
                return mc.normal_like(np.diff(np.log(np.maximum(NEARLY_ZERO, f[age_indices]))), 0.0, tau)
            priors += [smooth_rate]

        elif prior[0] == 'zero':
            tau_zero_rate = 1./(1e-4)**2
            
            age_start = int(prior[1])
            age_end = int(prior[2])
            age_indices = indices_for_range(age_mesh, age_start, age_end)
                               
            @mc.potential(name='zero_{%d,%d}^%s' % (age_start, age_end, rate))
            def zero_rate(f=rate, age_indices=age_indices, tau=tau_zero_rate):
                return mc.normal_like(f[age_indices], 0.0, tau)
            priors += [zero_rate]

        elif prior[0] == 'confidence':
            # prior only affects beta_binomial_rate model
            if not confidence_stoch:
                continue

            mu = float(prior[1])
            tau = float(prior[2])

            @mc.potential(name='prior_%s' % confidence_stoch)
            def confidence_potential(f=confidence_stoch, mu=mu, tau=tau):
                return mc.normal_like(f, mu, tau)
            priors += [confidence_potential]

        elif prior[0] == 'increasing':
            priors += derivative_sign_prior(rate, prior, deriv=1, sign=1)
        elif prior[0] == 'decreasing':
            priors += derivative_sign_prior(rate, prior, deriv=1, sign=-1)
        elif prior[0] == 'convex_down':
            priors += derivative_sign_prior(rate, prior, deriv=2, sign=-1)
        elif prior[0] == 'convex_up':
            priors += derivative_sign_prior(rate, prior, deriv=2, sign=1)

        elif prior[0] == 'unimodal':
            age_start = int(prior[1])
            age_end = int(prior[2])
            age_indices = indices_for_range(age_mesh, age_start, age_end)

            @mc.potential(name='unimodal_{%d,%d}^%s' % (age_start, age_end, rate))
            def unimodal_rate(f=rate, age_indices=age_indices, tau=1000.):
                df = np.diff(f[age_indices])
                sign_changes = pl.find((df[:-1] > NEARLY_ZERO) & (df[1:] < -NEARLY_ZERO))
                sign = np.ones(len(age_indices)-1)
                if len(sign_changes) > 0:
                    change_age = sign_changes[len(sign_changes)/2]
                    sign[change_age:] = -1.
                return -tau*np.dot(np.abs(df[:-1]), (sign * df[:-1] < 0))
            priors += [unimodal_rate]

        else:
            raise KeyError, 'Unrecognized prior: %s' % prior_str

    return priors

