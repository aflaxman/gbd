import numpy as np
import pymc as mc

import dismod3.settings
from dismod3.settings import MISSING, NEARLY_ZERO
from dismod3.utils import trim, clean, indices_for_range, rate_for_range, cscpm

import neg_binom_model as rate_model
import normal_model
import log_normal_model

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
      the observed data to use in the rate stoch likelihood functions
    
    Results
    -------
    vars : dict of PyMC stochs
      returns a dictionary of all the relevant PyMC objects for the
      generic disease model.
    """
    vars = {}


    param_type = 'all-cause_mortality'
    data = [d for d in data_list if d['data_type'] == 'all-cause mortality data']
    m_all_cause = dm.mortality(key % param_type, data)

    covariate_dict = dm.get_covariates()
    X_region, X_study = rate_model.regional_covariates(key, covariate_dict)
    est_mesh = dm.get_estimate_age_mesh()

    # update age_weights on non-incidence/prevalence data to reflect
    # prior prevalence distribution, if available
    prior_prev = dm.get_mcmc('emp_prior_mean', key % 'prevalence')
    if len(prior_prev) > 0:
        for d in data:
            if d['data_type'].startswith('incidence') or d['data_type'].startswith('prevalence'):
                continue
            age_indices = indices_for_range(est_mesh, d['age_start'], d['age_end'])
            d['age_weights'] = prior_prev[age_indices]
            d['age_weights'] /= sum(d['age_weights']) # age weights must sum to 1 (optimization of inner loop removed check on this)
                                      

    for param_type in ['incidence', 'remission', 'excess-mortality']:
        data = [d for d in data_list if d['data_type'] == '%s data' % param_type]

        lower_bound_data = []
        # TODO: include lower bound data when appropriate
        
        prior_dict = dm.get_empirical_prior(param_type)
        if prior_dict == {}:
            prior_dict.update(alpha=np.zeros(len(X_region)),
                              beta=np.zeros(len(X_study)),
                              gamma=-5*np.ones(len(est_mesh)),
                              sigma_alpha=[1.],
                              sigma_beta=[1.],
                              sigma_gamma=[1.],
                              # delta is filled in from the global prior dict in neg_binom setup
                              )
        vars[key % param_type] = rate_model.setup(dm, key % param_type, data,
                                                  emp_prior=prior_dict, lower_bound_data=lower_bound_data)

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
    def SC_0(C_0=C_0):
        return np.array([1. - C_0, C_0]).ravel()
    vars[key % 'bins'] = {'initial': [SC_0, C_0, logit_C_0]}
    
    
    # iterative solution to difference equations to obtain bin sizes for all ages
    import scipy.linalg
    import time;
    @mc.deterministic(name=key % 'bins')
    def SCpm(SC_0=SC_0, i=i, r=r, f=f, m_all_cause=m_all_cause, age_mesh=dm.get_param_age_mesh()):
        #t = time.time()
        SC = np.zeros([2, len(age_mesh)])
        p = np.zeros(len(age_mesh))
        m = np.zeros(len(age_mesh))
        
        SC[:,0] = SC_0
        p[0] = SC_0[1] / (SC_0[0] + SC_0[1])
        m[0] = trim(m_all_cause[age_mesh[0]] - f[age_mesh[0]] * p[0], .1*m_all_cause[age_mesh[0]], 1-NEARLY_ZERO)

        for ii, a in enumerate(age_mesh[:-1]):
            A = np.array([[-i[a]-m[ii],  r[a]          ],
                          [ i[a]     , -r[a]-m[ii]-f[a]]]) * (age_mesh[ii+1] - age_mesh[ii])

            SC[:,ii+1] = np.dot(scipy.linalg.expm(A), SC[:,ii])
            
            p[ii+1] = trim(SC[1,ii+1] / (SC[0,ii+1] + SC[1,ii+1]), NEARLY_ZERO, 1-NEARLY_ZERO)
            m[ii+1] = trim(m_all_cause[age_mesh[ii+1]] - f[age_mesh[ii+1]] * p[ii+1], .1*m_all_cause[age_mesh[ii+1]], 1-NEARLY_ZERO)

        SCpm = np.zeros([4, len(age_mesh)])
        SCpm[0:2,:] = SC
        SCpm[2,:] = p
        SCpm[3,:] = m
        return SCpm
        """
        """
        #ptime = time.time() - t
        #t = time.time()      
        cSCpm = cscpm(SC_0, i, r, f, m_all_cause, age_mesh, 1., NEARLY_ZERO)
        """
        ctime = time.time() - t
        print 'ptime =', ptime
        print 'ctime =', ctime
        print 'ptime/ctime =', ptime/ctime
        import pdb; pdb.set_trace()
        """
        return cSCpm

    vars[key % 'bins']['age > 0'] = [SCpm]

    
    # prevalence = # with condition / (# with condition + # without)
    @mc.deterministic(name=key % 'p')
    def p(SCpm=SCpm, param_mesh=dm.get_param_age_mesh(), est_mesh=dm.get_estimate_age_mesh()):
        return dismod3.utils.interpolate(param_mesh, SCpm[2,:], est_mesh)
    data = [d for d in data_list if d['data_type'] == 'prevalence data']
    prior_dict = dm.get_empirical_prior('prevalence')
    
    vars[key % 'prevalence'] = rate_model.setup(dm, key % 'prevalence', data, p, emp_prior=prior_dict)
    p = vars[key % 'prevalence']['rate_stoch']  # replace perfectly consistent p with version including level-bound priors
    
    # make a blank prior dict, to avoid weirdness
    blank_prior_dict = dict(alpha=np.zeros(len(X_region)),
                            beta=np.zeros(len(X_study)),
                            gamma=-5*np.ones(len(est_mesh)),
                            sigma_alpha=[1.],
                            sigma_beta=[1.],
                            sigma_gamma=[1.],
                            delta=100.,
                            sigma_delta=1.
                            )
    # cause-specific-mortality is a lower bound on p*f
    @mc.deterministic(name=key % 'pf')
    def pf(p=p, f=f):
        return (p+NEARLY_ZERO)*f
    # TODO: add a 'with-condition population mortality rate date' type
    # data = [d for d in data_list if d['data_type'] == 'with-condition population mortality rate data']
    data = []
    lower_bound_data = [d for d in data_list if d['data_type'] == 'cause-specific mortality data']
    vars[key % 'prevalence_x_excess-mortality'] = rate_model.setup(dm, key % 'pf', rate_stoch=pf, data_list=data, lower_bound_data=lower_bound_data, emp_prior=blank_prior_dict)
        

    # m = m_all_cause - f * p
    @mc.deterministic(name=key % 'm')
    def m(SCpm=SCpm, param_mesh=dm.get_param_age_mesh(), est_mesh=dm.get_estimate_age_mesh()):
        return dismod3.utils.interpolate(param_mesh,  SCpm[3,:], est_mesh)
    vars[key % 'm'] = m

    # m_with = m + f
    @mc.deterministic(name=key % 'm_with')
    def m_with(m=m, f=f):
        return m + f
    data = [d for d in data_list if d['data_type'] == 'mortality data']
    # TODO: test this
    #prior_dict = dm.get_empirical_prior('excess-mortality')  # TODO:  make separate prior for with-condition mortality
    vars[key % 'mortality'] = rate_model.setup(dm, key % 'm_with', data, m_with, emp_prior=blank_prior_dict)

    # mortality rate ratio = mortality with condition / mortality without
    @mc.deterministic(name=key % 'RR')
    def RR(m=m, m_with=m_with):
        return m_with / (m + .0001)
    data = [d for d in data_list if d['data_type'] == 'relative-risk data']
    vars[key % 'relative-risk'] = log_normal_model.setup(dm, key % 'relative-risk', data, RR)
    
    # standardized mortality rate ratio = mortality with condition / all-cause mortality
    @mc.deterministic(name=key % 'SMR')
    def SMR(m_with=m_with, m_all_cause=m_all_cause):
        return m_with / (m_all_cause + .0001)
    data = [d for d in data_list if d['data_type'] == 'smr data']
    vars[key % 'smr'] = log_normal_model.setup(dm, key % 'smr', data, SMR)

    # duration = E[time in bin C]
    @mc.deterministic(name=key % 'X')
    def X(r=r, m=m, f=f):
        hazard = r + m + f
        pr_not_exit = np.exp(-hazard)
        X = np.empty(len(hazard))
        X[-1] = 1 / hazard[-1]
        for i in reversed(range(len(X)-1)):
            X[i] = pr_not_exit[i] * (X[i+1] + 1) + 1 / hazard[i] * (1 - pr_not_exit[i]) - pr_not_exit[i]
        return X
    data = [d for d in data_list if d['data_type'] == 'duration data']
    vars[key % 'duration'] = normal_model.setup(dm, key % 'duration', data, X)

    # YLD[a] = disability weight * i[a] * X[a] * regional_population[a]
    @mc.deterministic(name=key % 'i*X')
    def iX(i=i, X=X, p=p, pop=rate_model.regional_population(key)):
        return i * X * (1-p) * pop 
    vars[key % 'incidence_x_duration'] = {'rate_stoch': iX}

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
            rate_model.store_mcmc_fit(dm, t, dm.vars[t])
