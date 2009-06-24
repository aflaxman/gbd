import numpy as np
import pymc as mc
import pymc.gp as gp

import dismod3.utils
from dismod3.utils import rate_for_range, generate_prior_potentials
from dismod3.settings import NEARLY_ZERO, MISSING

class LogitGPStep(mc.Metropolis):
    """ An attempt to speed up convergence of the MCMC for Beta Binomial models.

    Override the proposal method of the Metropolis StepMethod
    to draw from a realization of a carefully chosen Gaussian Process

    Parameters
    ----------
    rate_stoch : pymc.Stochastic
      a PyMC stochastic (or deterministic) object, with
      len(rate_stoch.value) == len(dm.get_estimation_age_mesh()).
      This is used to link beta-binomial stochs into a larger model,
      for example.

    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
      
    key : str
      the name of the key for everything about this model (priors,
      initial values, estimations)

    data_list : list of data dicts
      the observed data to use in the beta-binomial liklihood function

    verbose : optional
      Level of output verbosity: 0=none, 1=low, 2=medium, 3=high
    """
    def __init__(self, stochastic, dm, key, data_list, verbose=0):
        mc.Metropolis.__init__(self, stochastic, proposal_sd='LogitGP', verbose=verbose)
        
        M, C = dismod3.utils.uninformative_prior_gp(c=-5.,
                                                    diff_degree=2., amp=25., scale=200.)

        for d in data_list:
            if d['value'] == MISSING:
                print 'WARNING: data %d missing value' % id
                continue

            # ensure all rate data is valid
            a0 = d['age_start']
            a1 = d['age_end']
            a2 = .5 * (a0 + a1 + 1)
            assert a0 <= a1

            d_val = dm.value_per_1(d)
            if d_val < 0 or d_val > 1:
                print 'WARNING: data %d not in range (0,1)' % id
                continue
            elif d_val == 0:
                logit_val = -10.
            elif d_val == 1:
                logit_val = 10.
            else:
                logit_val = mc.logit(d_val)
        
            gp.observe(M, C,
                       [a for a in range(a0,a1)],
                       [logit_val for a in range(a0,a1)],
                       [100. for a in range(a0,a1)])

        prior_str = dm.get_priors(key)
        for p1 in prior_str.split(','):
            p = p1.split()
            if len(p) > 0 and p[0] == 'zero':
                a0 = int(p[1])
                a1 = int(p[2])
                gp.observe(M, C, [a for a in range(a0,a1+1)],
                           [-10. for a in range(a0,a1+1)],
                           [25. for a in range(a0,a1+1)])
        
        self.mesh = dm.get_param_age_mesh()
        self.M = M
        self.C = C

        if self.verbose:
            self.dm = dm
            self.data = data_list
            self.key = key

    def random(self):
        return gp.Realization(self.M, self.C)(self.mesh)
            
    def propose(self):
        """
        This method is called by step() to generate proposed values
        if self.proposal_distribution is "Normal" (i.e. no proposal specified).
        """
        random_gp = self.random()
        a = .125 * self.adaptive_scale_factor
        #a = .25

        if  mc.rbernoulli(.5):
            proposal = (1 - a) * self.stochastic.value + a * random_gp
        else:  # make chain reversible, simpler than computing Hasting Factor
            a = -a/(1-a)
            proposal = (1 - a) * self.stochastic.value + a * random_gp

        if self.verbose:
            try:
                import pylab as pl
                import time

                if self.key.find('incidence') != -1:
                    pl.subplot(1,3,2)
                elif self.key.find('case') != -1:
                    pl.subplot(1,3,1)
                else:
                    pl.subplot(1,3,3)
                

                dismod3.plotting.plot_intervals(self.dm, self.data)
                #dismod3.plotting.plot_prior(self.dm, self.key)
                #pl.clf()
                pl.plot(dismod3.utils.interpolate(
                        self.mesh,
                        mc.invlogit(random_gp),
                        range(100)),
                        alpha=.05, color='black', label='Random Sample', linewidth=3)

                pl.plot(dismod3.utils.interpolate(
                        self.mesh,
                        mc.invlogit(self.stochastic.value),
                        range(100)),
                        alpha=.05, color='blue', label='Current Value', linewidth=2)

                pl.plot(dismod3.utils.interpolate(
                        self.mesh,
                        mc.invlogit(proposal),
                        range(100)),
                        alpha=.05, color='red', label='Proposed Value', linewidth=2)
                #pl.legend()
                pl.figtext(.5,.5, 'a=%f' % a)
                pl.show()
                #time.sleep(1.)
            except ValueError:
                pass

        self.stochastic.value = proposal
    

    def hastings_factor(self):
        return 0. #self.gp_logp_p - self.gp_logp
