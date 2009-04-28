import numpy as np
import pymc as mc

from bayesian_models import probabilistic_utils
import beta_binomial_model as rate_model

output_data_types = ['incidence data', 'remission data', 'case fatality data', 'prevalence data', 'duration data']

def initialize(dm):
    for data_type in ['incidence data', 'remission data', 'case fatality data']:
        rate_model.initialize(dm, data_type)
    dm.set_units('prevalence data', '(per person)')
    dm.set_units('duration data', '(years)')

def map_fit(dm, vars):
    map = mc.MAP([v.values() for v in vars.values()])
    map.fit(method='fmin_powell', iterlim=500, tol=.001, verbose=1)
    for data_type in output_data_types:
        dm.set_map(data_type, vars[data_type]['rate_stoch'].value)
    return map

def mcmc_fit(dm, vars, data_type='prevalence data'):
    mcmc = mc.MCMC([v.values() for v in vars.values()])
    mcmc.sample(iter=400000, burn=100000, thin=200, verbose=1)
    for data_type in output_data_types:
        rate_model.store_mcmc_fit(dm, vars[data_type]['rate_stoch'], data_type)
    

def setup(dm):
    """
    Generate the PyMC variables for the generic disease
    model, and return it as a dictionary of vars and lists of vars
    """
    vars = {}
    
    vars['incidence data'] = rate_model.setup(dm, 'incidence data')
    i = vars['incidence data']['rate_stoch']

    vars['remission data'] = rate_model.setup(dm, 'remission data')
    r = vars['remission data']['rate_stoch']

    vars['case fatality data'] = rate_model.setup(dm, 'case fatality data')
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
    vars['prevalence data'] = rate_model.setup(dm, 'prevalence data', p)
    
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
