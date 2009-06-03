import numpy as np
import pymc as mc

from bayesian_models import probabilistic_utils
from dismod3.utils import clean

import beta_binomial_model as rate_model

output_data_types = ['incidence', 'remission', 'case-fatality', 'prevalence', 'duration']

def fit(dm, method='map'):
    """ Generate an estimate of the generic disease model parameters
    using maximum a posteriori liklihood (MAP) or Markov-chain Monte
    Carlo (MCMC)

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)

    method : string, optional
      the parameter estimation method, either 'map' or 'mcmc'

    Example
    -------
    >>> import dismod3
    >>> import dismod3.generic_disease_model as model
    >>> dm = dismod3.get_disease_model(1)
    >>> model.fit(dm, method='map')
    >>> model.fit(dm, method='mcmc')
    """
    if not hasattr(dm, 'vars'):
        for param_type in ['incidence', 'remission', 'case-fatality']:
            # find initial values for these rates
            data =  [d for d in dm.data if clean(d['data_type']).find(param_type) != -1]

            # use a random subset of the data if there is a lot of it,
            # to speed things up
            if len(data) > 25:
                dm.fit_initial_estimate(param_type, random.sample(data,25))
            else:
                dm.fit_initial_estimate(param_type, data)

            dm.set_units(param_type, '(per person-year)')

        dm.set_units('prevalence', '(per person)')
        dm.set_units('duration', '(years)')

        dm.vars = setup(dm, dm.data)

    if method == 'map':
        if not hasattr(dm, 'map'):
            dm.map = mc.MAP(dm.vars)
            
        dm.map.fit(method='fmin_powell', iterlim=500, tol=.001, verbose=1)
        for t in output_data_types:
            dm.set_map(t, dm.vars[t]['rate_stoch'].value)
    elif method == 'mcmc':
        if not hasattr(dm, 'mcmc'):
            dm.mcmc = mc.MCMC(dm.vars)
            for key in dm.vars:
                stochs = dm.vars[key].get('logit_p_stochs', [])
                if len(stochs) > 0:
                    dm.mcmc.use_step_method(mc.AdaptiveMetropolis, stochs)

        dm.mcmc.sample(iter=60*1000, burn=10*1000, thin=50, verbose=1)
        for t in output_data_types:
            rate_model.store_mcmc_fit(dm, t, dm.vars[t]['rate_stoch'])


def setup(dm, key='%s', data_list=None):
    """ Generate the PyMC variables for a generic disease model

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)

    key : str, optional
      a string for modifying the names of the stochs in this model,
      must contain a single %s that will be substituted

    data_list : list of data dicts
      the observed data to use in the beta-binomial liklihood function
    
    Results
    -------
    vars : dict of PyMC stochs
      returns a dictionary of all the relevant PyMC objects for the
      generic disease model.
    """

    if data_list == None:
        data_list = dm.data
    
    vars = {}
    
    for param_type in ['incidence', 'remission', 'case-fatality']:
        # find initial values for these rates
        data = [d for d in data_list if clean(d['data_type']).find(param_type) != -1]
        vars[key % param_type] = rate_model.setup(dm, key % param_type, data)

    i = vars[key % 'incidence']['rate_stoch']
    r = vars[key % 'remission']['rate_stoch']
    f = vars[key % 'case-fatality']['rate_stoch']

    param_type = 'all-cause_mortality'
    data = [d for d in data_list if clean(d['data_type']).find(param_type) != -1]

    m = dm.mortality(key % param_type, data)
    
    # TODO: make error in C_0 a semi-informative stochastic variable
    logit_C_0 = mc.Normal('logit(C_0^%s)' % key, -5., 1.e-2)
    @mc.deterministic
    def C_0(logit_C_0=logit_C_0):
        return mc.invlogit(logit_C_0)
    
    @mc.deterministic
    def S_0(C_0=C_0):
        return max(0.0, 1.0 - C_0)
    vars[key % 'bins'] = {'initial': [S_0, C_0, logit_C_0]}
    
    # iterative solution to difference equations to obtain bin sizes for all ages
    age_len = len(dm.get_estimate_age_mesh())
    @mc.deterministic
    def S_C_D_M(S_0=S_0, C_0=C_0, i=i, r=r, f=f, m=m, age_len=age_len):
        S = np.zeros(age_len)
        C = np.zeros(age_len)
        D = np.zeros(age_len)
        M = np.zeros(age_len)
        
        S[0] = S_0
        C[0] = C_0
        D[0] = 0.0
        M[0] = 0.0
        
        for a in range(age_len - 1):
            S[a+1] = S[a]*max(0, 1-i[a]-m[a]) + C[a]*r[a]
            C[a+1] = S[a]*i[a]                + C[a]*max(0, 1-r[a]-m[a]-f[a])
            D[a+1] =                            C[a]*f[a]                     + D[a]
            M[a+1] = S[a]*m[a]                + C[a]*m[a]                            + M[a]
                
        return S,C,D,M
    vars[key % 'bins']['age > 0'] = [S_C_D_M]

    # prevalence = # with condition / (# with condition + # without)
    @mc.deterministic
    def p(S_C_D_M=S_C_D_M, tau_p=1./.01**2):
        S,C,D,M = S_C_D_M
        return probabilistic_utils.trim(C/(S+C),
                                        probabilistic_utils.NEARLY_ZERO,
                                        1. - probabilistic_utils.NEARLY_ZERO)
    
    data = [d for d in data_list if clean(d['data_type']).find('prevalence') != -1]
    vars[key % 'prevalence'] = rate_model.setup(dm, key % 'prevalence', data, p)
    
    # duration = E[time in bin C]
    @mc.deterministic
    def X(r=r, m=m, f=f):
        pr_exit = 1 - r - m - f
        X = np.empty(len(pr_exit))
        t = 1.0
        for i in xrange(len(X)-1,-1,-1):
            X[i] = t*pr_exit[i]
            t = 1+X[i]
        return X
    vars[key % 'duration'] = {'rate_stoch': X}

    return vars
