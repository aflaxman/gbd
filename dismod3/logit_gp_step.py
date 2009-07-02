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
        
        M, C = dismod3.utils.uninformative_prior_gp(c=0., diff_degree=2., amp=25., scale=200.)

        min_val = min([1.e-2] + [dm.value_per_1(d) for d in data_list if dm.value_per_1(d) > 0])
        
        for d in data_list:
            try:
                age_indices, age_weights, logit_val, logit_se = values_from(dm, d, min_val)
            except ValueError:
                continue

            #TODO:  combine overlapping age intervals more systematically
            gp.observe(M, C,
                       [a for a in age_indices],
                       [logit_val for a in age_indices],
                       [25. for a in age_indices])

        prior_str = dm.get_priors(key)
        for p1 in prior_str.split(','):
            p = p1.split()
            if len(p) > 0 and p[0] == 'zero':
                a0 = int(p[1])
                a1 = int(p[2])
                gp.observe(M, C, [a for a in range(a0,a1+1)],
                           [mc.logit(min_val / 100) for a in range(a0,a1+1)],
                           [25. for a in range(a0,a1+1)])

        self.scale = .0125
        
        self.mesh = dm.get_param_age_mesh()
        self.M = M
        self.C = C

    def random(self):
        return gp.Realization(self.M, self.C)(self.mesh)
            
    def propose(self):
        """
        This method is called by step() to generate proposed values
        if self.proposal_distribution is "Normal" (i.e. no proposal specified).
        """
        random_gp = self.random()

        if self.scale < 1.:
            a = self.scale * self.adaptive_scale_factor
            if  mc.rbernoulli(.5):
                proposal = (1 - a) * self.stochastic.value + a * random_gp
            else:  # make chain reversible, simpler than computing Hasting Factor
                a = -a/(1-a)
                proposal = (1 - a) * self.stochastic.value + a * random_gp
        else: # hacky way to draw purely from random gp
            proposal = random_gp

        self.stochastic.value = proposal
