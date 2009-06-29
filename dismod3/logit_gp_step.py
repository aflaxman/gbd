import numpy as np
import pymc as mc
import pymc.gp as gp

import dismod3.utils
from dismod3.utils import debug, rate_for_range, generate_prior_potentials
from dismod3.settings import NEARLY_ZERO, MISSING
from dismod3.logit_normal_model import values_from

class LogitGPStep(mc.Metropolis):
    """ An attempt to speed up convergence of the MCMC for "random effect logistic models".

    Override the proposal method of the Metropolis StepMethod
    to draw from a realization of a carefully chosen Gaussian Process

    Parameters
    ----------
    stochastic : pymc.Stochastic

    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
      
    key : str
      the name of the key for everything about this model (priors,
      initial values, estimations)

    data_list : list of data dicts
      the observed data to use in the approximation of the prior distribution

    verbose : optional
      Level of output verbosity: 0=none, 1=low, 2=medium, 3=high
    """
    def __init__(self, stochastic, dm, key, data_list, verbose=0):
        mc.Metropolis.__init__(self, stochastic, proposal_sd='LogitGP', verbose=verbose)
        
        M, C = dismod3.utils.uninformative_prior_gp(c=5.,
                                                    diff_degree=2., amp=25., scale=200.)

        for d in data_list:
            try:
                age_indices, age_weights, logit_val, logit_se = values_from(dm, d)
            except ValueError:
                continue

            gp.observe(M, C,
                       [a for a in age_indices],
                       [logit_val for a in age_indices],
                       [100. for a in age_indices])

        prior_str = dm.get_priors(key)
        for p1 in prior_str.split(','):
            p = p1.split()
            if len(p) > 0 and p[0] == 'zero':
                a0 = int(p[1])
                a1 = int(p[2])
                gp.observe(M, C, [a for a in range(a0,a1+1)],
                           [-10. for a in range(a0,a1+1)],
                           [25. for a in range(a0,a1+1)])

        self.scale = .125
        
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
        a = self.scale * self.adaptive_scale_factor

        if a < 1.:
            if  mc.rbernoulli(.5):
                proposal = (1 - a) * self.stochastic.value + a * random_gp
            else:  # make chain reversible, simpler than computing Hasting Factor
                a = -a/(1-a)
                proposal = (1 - a) * self.stochastic.value + a * random_gp
        else: # hacky way to draw purely from random gp
            proposal = random_gp
            
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
