import numpy as np
import pymc as mc

import dismod3.settings
from dismod3.settings import MISSING, NEARLY_ZERO
from dismod3.utils import trim, clean, indices_for_range, rate_for_range

import neg_binom_model as rate_model
import normal_model

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
      the observed data to use in the rate stoch likelihood functions

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

    m_all_cause = dm.mortality(key % param_type, data)
    
    for param_type in ['incidence', 'remission', 'excess-mortality']:
        data = [d for d in data_list if clean(d['data_type']).find(param_type) != -1]
        prior_dict = dm.get_empirical_prior(param_type)
        vars[key % param_type] = rate_model.setup(dm, key % param_type, data, emp_prior=prior_dict)

    i = vars[key % 'incidence']['rate_stoch']
    r = vars[key % 'remission']['rate_stoch']
    f = vars[key % 'excess-mortality']['rate_stoch']

    # Initial population with condition
    logit_C_0 = mc.Normal('logit_%s' % (key % 'C_0'), -5., 1., value=-5.)
    @mc.deterministic(name=key % 'C_0')
    def C_0(logit_C_0=logit_C_0):
        return mc.invlogit(logit_C_0)
    
    # Initial population without condition
    @mc.deterministic(name=key % 'S_0')
    def S_0(C_0=C_0):
        return 1. - C_0
    vars[key % 'bins'] = {'initial': [S_0, C_0, logit_C_0]}
    
    # iterative solution to difference equations to obtain bin sizes for all ages
    age_len = len(dm.get_estimate_age_mesh())
    @mc.deterministic(name=key % 'bins')
    def S_C_D_M_p_m(S_0=S_0, C_0=C_0, i=i, r=r, f=f, m_all_cause=m_all_cause, age_len=age_len):
        import scipy.linalg
        
        SCDM = np.zeros([4, age_len])
        p = np.zeros(age_len)
        m = np.zeros(age_len)

        SCDM[0,0] = S_0
        SCDM[1,0] = C_0
        SCDM[2,0] = NEARLY_ZERO
        SCDM[3,0] = NEARLY_ZERO

        p[0] = SCDM[1,0] / (SCDM[0,0] + SCDM[1,0] + NEARLY_ZERO)
        m[0] = trim(m_all_cause[0] - f[0] * p[0], NEARLY_ZERO, 1-NEARLY_ZERO)
        
        for a in range(age_len - 1):
            A = [[-i[a]-m[a],  r[a]          , 0., 0.],
                 [ i[a]     , -r[a]-m[a]-f[a], 0., 0.],
                 [      m[a],       m[a]     , 0., 0.],
                 [        0.,            f[a], 0., 0.]]

#             SCDM[:,a+1] = np.dot(np.eye(4) + A, SCDM[:,a])
#             SCDM[:,a+1] = np.dot(np.eye(4) + A + .5*np.dot(A,A), SCDM[:,a])
            SCDM[:,a+1] = np.dot(scipy.linalg.expm(A), SCDM[:,a])
#             SCDM[:,a+1] = np.dot(scipy.linalg.expm3(A,2), SCDM[:,a])
            
            p[a+1] = SCDM[1,a+1] / (SCDM[0,a+1] + SCDM[1,a+1] + NEARLY_ZERO)
            m[a+1] = trim(m_all_cause[a+1] - f[a+1] * p[a+1], .1*m_all_cause[a+1], 1-NEARLY_ZERO)

        SCDMpm = np.zeros([6, age_len])
        SCDMpm[0:4,:] = SCDM
        SCDMpm[4,:] = np.maximum(p, NEARLY_ZERO)
        SCDMpm[5,:] = np.maximum(m, NEARLY_ZERO)
        
        return SCDMpm
    vars[key % 'bins']['age > 0'] = [S_C_D_M_p_m]

    # prevalence = # with condition / (# with condition + # without)
    @mc.deterministic(name=key % 'p')
    def p(SCDMpm=S_C_D_M_p_m):
        return SCDMpm[4,:]
    data = [d for d in data_list if clean(d['data_type']).find('prevalence') != -1]
    prior_dict = dm.get_empirical_prior('prevalence')

    vars[key % 'prevalence'] = rate_model.setup(dm, key % 'prevalence', data, p, emp_prior=prior_dict)
    
    # m = m_all_cause - f * p
    @mc.deterministic(name=key % 'm')
    def m(SCDMpm=S_C_D_M_p_m):
        return SCDMpm[5,:]
    vars[key % 'm'] = m

    # m_with = m + f
    @mc.deterministic(name=key % 'm_with')
    def m_with(m=m, f=f):
        return m + f
    data = [d for d in data_list if clean(d['data_type']).find('mortality') != -1]
    prior_dict = dm.get_empirical_prior('excess-mortality')  # TODO:  make separate prior for with-condition mortality
    vars[key % 'mortality'] = rate_model.setup(dm, key % 'm_with', data, m_with, emp_prior=prior_dict)

    # relative risk = mortality with condition / mortality without
    @mc.deterministic(name=key % 'RR')
    def RR(m=m, f=f):
        return (m + f) / m
    data = [d for d in data_list if clean(d['data_type']).find('relative-risk') != -1]
    vars[key % 'relative-risk'] = normal_model.setup(dm, key % 'relative-risk', data, RR)
    
    # duration = E[time in bin C]
    @mc.deterministic(name=key % 'X')
    def X(r=r, m=m, f=f):
        pr_exit = np.exp(- r - m - f)
        X = np.empty(len(pr_exit))
        t = 1.
        for i in xrange(len(X)-1,-1,-1):
            X[i] = t*pr_exit[i]
            t = 1+X[i]
        return X
    data = [d for d in data_list if clean(d['data_type']).find('duration') != -1]
    vars[key % 'duration'] = normal_model.setup(dm, key % 'duration', data, X)

    # YLD[a] = disability weight * i[a] * X[a] * regional_population[a]
    @mc.deterministic(name=key % 'i*X')
    def iX(i=i, X=X):
        return i * X
    vars[key % 'incidence_x_duration'] = {'rate_stoch': iX}

    # SMR

    return vars









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
        for param_type in ['incidence', 'remission', 'excess-mortality']:
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
