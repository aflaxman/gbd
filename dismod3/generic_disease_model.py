import numpy as np
import pymc as mc

import dismod3.settings
from dismod3.settings import MISSING, NEARLY_ZERO
from dismod3.utils import trim, clean, indices_for_range, rate_for_range


## alternative rate models  (pick one)
import logit_normal_model as rate_model
#import beta_binomial_model as rate_model

import normal_model

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

        dm.vars = setup(dm)

    if method == 'map':
        if not hasattr(dm, 'map'):
            dm.map = mc.MAP(dm.vars)
            
        try:
            dm.map.fit(method='fmin_powell', iterlim=500, tol=.001, verbose=1)
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass

        for t in dismod3.settings.output_data_types:
            t = clean(t)
            val = dm.vars[t]['rate_stoch'].value
            dm.set_map(t, val)
            dm.set_initial_value(t, val)  # better initial value may save time in the future
    elif method == 'mcmc':
        if not hasattr(dm, 'mcmc'):
            dm.mcmc = mc.MCMC(dm.vars)
            for key in dm.vars:
                stochs = dm.vars[key].get('logit_p_stochs', [])
                if len(stochs) > 0:
                    dm.mcmc.use_step_method(mc.AdaptiveMetropolis, stochs)

        try:
            dm.mcmc.sample(iter=60*1000, burn=10*1000, thin=50, verbose=1)
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass
        for t in dismod3.settings.output_data_types:
            t = clean(t)
            rate_model.store_mcmc_fit(dm, t, dm.vars[t]['rate_stoch'])

def setup(dm, key='%s', data_list=None, regional_population=None):
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

    regional_population : list, optional
      the population of the region, on dm.get_estimate_age_mesh(), for calculating YLDs by age
    
    Results
    -------
    vars : dict of PyMC stochs
      returns a dictionary of all the relevant PyMC objects for the
      generic disease model.
    """

    if data_list == None:
        data_list = dm.data
    
    vars = {}

    param_type = 'all-cause_mortality'
    data = [d for d in data_list if clean(d['data_type']).find(param_type) != -1]
    m = dm.mortality(key % param_type, data)
    
    for param_type in ['incidence', 'remission', 'case-fatality']:
        data = [d for d in data_list if clean(d['data_type']).find(param_type) != -1]

        if param_type == 'case-fatality':
            # apply relative-risk data to case-fatality rate directly
            est_mesh = dm.get_estimate_age_mesh()
            rr_data = [d for d in data_list if clean(d['data_type']).find('relative-risk') != -1]
            for d in rr_data:
                d_val = dm.value_per_1(d)
                if d_val != MISSING:
                    age_indices = indices_for_range(est_mesh, d['age_start'], d['age_end'])
                    age_weights = d['age_weights']
                    
                    d['data_type'] = 'case-fatality data'
                    # case-fatality = all-cause mortality * (relative risk - 1)
                    d['value'] = (d['value'] - 1) * rate_for_range(m, age_indices, age_weights)
                    d['units'] = 'per 1.0'

                    # TODO: calculate standard error of cf from rr
                    d['standard_error'] = MISSING

                    data.append(d)
                    
        prior_dict = dm.get_empirical_prior(param_type)
        vars[key % param_type] = rate_model.setup(dm, key % param_type, data, emp_prior=prior_dict)

    i = vars[key % 'incidence']['rate_stoch']
    r = vars[key % 'remission']['rate_stoch']
    f = vars[key % 'case-fatality']['rate_stoch']

    # relative risk = mortality with condition / mortality without
    @mc.deterministic(name='RR_%s' % key)
    def RR(m=m, f=f):
        return (m + f) / m
    # convert relative risk data into case-fatality data, and use it there, but permit
    # priors on relative risk curve directly
    # data = [d for d in data_list if clean(d['data_type']).find('relative-risk') != -1]
    data = []
    vars[key % 'relative-risk'] = normal_model.setup(dm, key % 'relative-risk', data, RR)
    
    # TODO: make error in C_0 a semi-informative stochastic variable
    logit_C_0 = mc.Normal('logit(%s)' % (key % 'C_0'), -5., 1.e-2)
    @mc.deterministic(name=key % 'C_0')
    def C_0(logit_C_0=logit_C_0):
        return mc.invlogit(logit_C_0)
    
    @mc.deterministic(name=key % 'S_0')
    def S_0(C_0=C_0):
        return 1. - C_0
    vars[key % 'bins'] = {'initial': [S_0, C_0, logit_C_0]}
    
    # iterative solution to difference equations to obtain bin sizes for all ages
    age_len = len(dm.get_estimate_age_mesh())
    @mc.deterministic(name='bins_%s' % key)
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
    @mc.deterministic(name='p_%s' % key)
    def p(S_C_D_M=S_C_D_M, tau_p=1./.01**2):
        S,C,D,M = S_C_D_M
        return trim(C / (S + C + NEARLY_ZERO), NEARLY_ZERO, 1. - NEARLY_ZERO)
    data = [d for d in data_list if clean(d['data_type']).find('prevalence') != -1]

    prior_dict = dm.get_empirical_prior('prevalence')
    vars[key % 'prevalence'] = rate_model.setup(dm, key % 'prevalence', data, p, emp_prior=prior_dict)
    
    # duration = E[time in bin C]
    @mc.deterministic(name='X_%s' % key)
    def X(r=r, m=m, f=f):
        pr_exit = 1 - r - m - f
        X = np.empty(len(pr_exit))
        t = 1.0
        for i in xrange(len(X)-1,-1,-1):
            X[i] = t*pr_exit[i]
            t = 1+X[i]
        return X
    vars[key % 'duration'] = {'rate_stoch': X}

    # YLD[a] = i[a] * X[a] * regional_population[a]
    if regional_population == None:
        regional_population = np.ones(len(dm.get_estimate_age_mesh()))
    @mc.deterministic(name='YLD_%s' % key)
    def YLD(i=i, X=X, pop=regional_population):
        return i * X * pop
    vars[key % 'yld'] = {'rate_stoch': YLD}

    return vars
