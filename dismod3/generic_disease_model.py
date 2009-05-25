import numpy as np
import pymc as mc

from bayesian_models import probabilistic_utils
import beta_binomial_model as rate_model

output_data_types = ['incidence data', 'remission data', 'case fatality data', 'prevalence data', 'duration data']

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
    >>> dm = dismod3.get_disease_model(854)
    >>> model.fit(dm, method='map')
    >>> model.fit(dm, method='mcmc')
    """
    if not hasattr(dm, 'vars'):
        initialize(dm)

    if method == 'map':
        if not hasattr(dm, 'map'):
            dm.map = mc.MAP(dm.vars)
        dm.map.fit(method='fmin_powell', iterlim=500, tol=.001, verbose=1)
        for t in output_data_types:
            dm.set_map(t, dm.vars[t]['rate_stoch'].value)
    elif method == 'mcmc':
        if not hasattr(dm, 'mcmc'):
            dm.mcmc = mc.MCMC(dm.vars)
            for est_vars in dm.vars.values():
                if len(est_vars.get('logit_p_stochs', [])) > 0:
                    dm.mcmc.use_step_method(
                        mc.AdaptiveMetropolis, est_vars['logit_p_stochs'])
                    
        dm.mcmc.sample(iter=60*1000, burn=10*1000, thin=50, verbose=1)
        for t in output_data_types:
            rate_model.store_mcmc_fit(dm, dm.vars[t]['rate_stoch'], t)



def initialize(dm):
    """ Initialize the stochastic and deterministic random variables
    for the generic disease model

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
    
    Results
    -------
    * Sets dm's param_age_mesh and estimate_age_mesh, if they are not
      already set.

    * Sets the units of all estimates in the dm

    * Create PyMC variables for the generic disease model, and store
      them in dm.vars
    """
    if dm.get_param_age_mesh() == []:
        dm.set_param_age_mesh([0.0, 10.0, 20.0, 30.0, 40.0,
                               50.0, 60.0, 70.0, 80.0, 90.0, 100.0])
    if dm.get_estimate_age_mesh() == []:
        dm.set_estimate_age_mesh(range(MAX_AGE))

    # sort the data by data types
    dm.data_by_type = {}
    for d in dm.data:
        t = d['data_type']
        dm.data_by_type[t] = dm.data_by_type.get(t, []) + [d]

    # find initial values for the rates that can be set
    for data_type in ['incidence data', 'remission data', 'case fatality data']:
        rate_model.initialize(dm, data_type)
    dm.set_units('prevalence data', '(per person)')
    dm.set_units('duration data', '(years)')

    dm.vars = setup(dm, data_type)


def setup(dm, data_type='prevalence data'):
    """ Generate the PyMC variables for a generic disease
    model

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
    
    Results
    -------
    vars : dict of PyMC stochs
      returns a dictionary of all the relevant PyMC objects for the
      generic disease model.
    """
    
    vars = {}
    
    vars['incidence data'] = rate_model.setup(
        dm, dm.data_by_type.get('incidence data', []), 'incidence data')
    i = vars['incidence data']['rate_stoch']

    vars['remission data'] = rate_model.setup(
        dm, dm.data_by_type.get('remission data', []), 'remission data')
    r = vars['remission data']['rate_stoch']

    vars['case fatality data'] = rate_model.setup(
        dm, dm.data_by_type.get('case fatality data',[]), 'case fatality data')
    f = vars['case fatality data']['rate_stoch']

    m = dm.mortality()
    
    # TODO: make error in C_0 a semi-informative stochastic variable
    logit_C_0 = mc.Normal('logit(C_0)', 0., 1.e-2)
    @mc.deterministic
    def C_0(logit_C_0=logit_C_0):
        return mc.invlogit(logit_C_0)
    
    @mc.deterministic
    def S_0(C_0=C_0):
        return max(0.0, 1.0 - C_0)
    vars['bins'] = {'initial': [S_0, C_0, logit_C_0]}
    
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
            S[a+1] = S[a]*(1-i[a]-m[a]) + C[a]*r[a]
            C[a+1] = S[a]*i[a]          + C[a]*(1-r[a]-m[a]-f[a])
            D[a+1] =                      C[a]*f[a]               + D[a]
            M[a+1] = S[a]*m[a]          + C[a]*m[a]                      + M[a]
                
        return S,C,D,M
    vars['bins']['age > 0'] = [S_C_D_M]

    # prevalence = # with condition / (# with condition + # without)
    @mc.deterministic
    def p(S_C_D_M=S_C_D_M, tau_p=1./.01**2):
        S,C,D,M = S_C_D_M
        return probabilistic_utils.trim(C/(S+C),
                                        probabilistic_utils.NEARLY_ZERO,
                                        1. - probabilistic_utils.NEARLY_ZERO)
    vars['prevalence data'] = rate_model.setup(dm, dm.data_by_type['prevalence data'],
                                               'prevalence data', p)
    
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
    vars['duration data'] = {'rate_stoch': X}

    return vars
