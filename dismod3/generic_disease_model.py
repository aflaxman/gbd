import pylab as pl
import pymc as mc

import dismod3

import neg_binom_model
import normal_model
import log_normal_model

def setup(dm, key='%s+north_america_high_income+2005+male', data_list=None):
    """ Generate the PyMC variables for a generic disease model

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)

    key : str, optional
      a string for modifying the names of the stochs in this model,
      must contain a single %s that will be substituted

    data_list : list of data dicts, optional
      the observed data to use in the rate stoch likelihood functions
    
    Results
    -------
    vars : dict of PyMC stochs
      returns a dictionary of all the relevant PyMC objects for the
      generic disease model.
    """
    if not data_list:
        data_list = dm.data

    type, region, sex, year = dismod3.utils.type_region_year_sex_from_key(key)

    vars = {}

    # setup all-cause mortality 
    param_type = 'all-cause_mortality'
    m_all_cause = dm.mortality(key % param_type)
    vars[key%param_type] = m_all_cause

    # make covariate vectors and estimation vectors to know dimensions of these objects
    covariate_dict = dm.get_covariates()
    derived_covariate = dm.get_derived_covariate_values()
    X_region, X_study = neg_binom_model.regional_covariates(key, covariate_dict, derived_covariate)
    est_mesh = dm.get_estimate_age_mesh()

    # update age_weights on non-incidence/prevalence data to reflect
    # prior prevalence distribution, if available
    prior_prev = dm.get_mcmc('emp_prior_mean', key % 'prevalence')
    if len(prior_prev) > 0:
        for d in data_list:
            if d['data_type'].startswith('incidence') or d['data_type'].startswith('prevalence'):
                continue
            age_indices = dismod3.utils.indices_for_range(est_mesh, d['age_start'], d['age_end'])
            d['age_weights'] = prior_prev[age_indices]
            d['age_weights'] /= sum(d['age_weights']) # age weights must sum to 1 (optimization of inner loop removed check on this)
                                      
    # create negative binomial models for incidence, remission, and
    # excess-mortality (which are all treated as "free" parameters)
    for param_type in ['incidence', 'remission', 'excess-mortality']:
        data = [d for d in data_list if d['data_type'] == '%s data' % param_type]

        lower_bound_data = [] # TODO: include lower bound data when appropriate (this has not come up yet)
        
        prior_dict = dm.get_empirical_prior(param_type)  # use empirical priors for the type/region/year/sex if available
        if prior_dict == {} and region != 'world':  # otherwise use weakly informative priors
            prior_dict.update(alpha=pl.zeros(len(X_region)),
                              beta=pl.zeros(len(X_study)),
                              gamma=-5*pl.ones(len(est_mesh)),
                              sigma_alpha=[1.],
                              sigma_beta=[1.],
                              sigma_gamma=[10.],
                              # delta is filled in from the global prior dict in neg_binom setup
                              )
        vars[key % param_type] = neg_binom_model.setup(dm, key % param_type, data,
                                                  emp_prior=prior_dict, lower_bound_data=lower_bound_data)

    # create nicer names for the rate stochastic from each neg-binom rate model
    i = vars[key % 'incidence']['rate_stoch']
    r = vars[key % 'remission']['rate_stoch']
    f = vars[key % 'excess-mortality']['rate_stoch']

    # initial fraction of population with the condition
    logit_C_0 = mc.Normal('logit_%s' % (key % 'C_0'), -5., 10.**-2, value=-5.)  # represet C_0 in logit space to allow unconstrained posterior maximization
    @mc.deterministic(name=key % 'C_0')
    def C_0(logit_C_0=logit_C_0):
        return mc.invlogit(logit_C_0)
    
    # initial fraction population with and without condition
    @mc.deterministic(name=key % 'SC_0')
    def SC_0(C_0=C_0):
        return pl.array([1. - C_0, C_0]).ravel()
    vars[key % 'bins'] = {'initial': dict(SC_0=SC_0, C_0=C_0, logit_C_0=logit_C_0)}
    
    
    # iterative solution to difference equations to obtain bin sizes for all ages
    import scipy.linalg
    @mc.deterministic(name=key % 'bins')
    def SCpm(SC_0=SC_0, i=i, r=r, f=f, m_all_cause=m_all_cause, age_mesh=dm.get_param_age_mesh()):
        SC = pl.zeros([2, len(age_mesh)])
        p = pl.zeros(len(age_mesh))
        m = pl.zeros(len(age_mesh))
        
        SC[:,0] = SC_0
        p[0] = SC_0[1] / (SC_0[0] + SC_0[1])
        m[0] = dismod3.utils.trim(m_all_cause[age_mesh[0]] - f[age_mesh[0]] * p[0],
                                  .1*m_all_cause[age_mesh[0]],
                                  1-dismod3.settings.NEARLY_ZERO)  # trim m[0] to avoid numerical instability

        for ii, a in enumerate(age_mesh[:-1]):
            A = pl.array([[-i[a]-m[ii], r[a]           ],
                          [ i[a]     , -r[a]-m[ii]-f[a]]]) * (age_mesh[ii+1] - age_mesh[ii])

            SC[:,ii+1] = pl.dot(scipy.linalg.expm(A), SC[:,ii])
            
            p[ii+1] = dismod3.utils.trim(SC[1,ii+1] / (SC[0,ii+1] + SC[1,ii+1]),
                                         dismod3.settings.NEARLY_ZERO,
                                         1-dismod3.settings.NEARLY_ZERO)
            m[ii+1] = dismod3.utils.trim(m_all_cause[age_mesh[ii+1]] - f[age_mesh[ii+1]] * p[ii+1],
                                         .1*m_all_cause[age_mesh[ii+1]],
                                         pl.inf)

        SCpm = pl.zeros([4, len(age_mesh)])
        SCpm[0:2,:] = SC
        SCpm[2,:] = p
        SCpm[3,:] = m
        return SCpm

    vars[key % 'bins']['age > 0'] = [SCpm]

    
    # prevalence = # with condition / (# with condition + # without)
    @mc.deterministic(name=key % 'p')
    def p(SCpm=SCpm, param_mesh=dm.get_param_age_mesh(), est_mesh=dm.get_estimate_age_mesh()):
        return dismod3.utils.interpolate(param_mesh, SCpm[2,:], est_mesh, kind='linear')
    data = [d for d in data_list if d['data_type'] == 'prevalence data']
    prior_dict = dm.get_empirical_prior('prevalence')
    # Does it ever help to over-ride the fully bayesian prior with this alternative expert prior?
    # if prior_dict == {} and region != 'world':
    #     prior_dict.update(alpha=pl.zeros(len(X_region)),
    #                       beta=pl.zeros(len(X_study)),
    #                       gamma=-5*pl.ones(len(est_mesh)),
    #                       sigma_alpha=[1.],
    #                       sigma_beta=[1.],
    #                       sigma_gamma=[10.],
    #                       # delta is filled in from the global prior dict in neg_binom setup
    #                       )
    
    vars[key % 'prevalence'] = neg_binom_model.setup(dm, key % 'prevalence', data, p, emp_prior=prior_dict)
    p = vars[key % 'prevalence']['rate_stoch']  # replace perfectly consistent p with version including level-bound priors
    
    # make a blank prior dict, to avoid weirdness
    blank_prior_dict = dict(alpha=pl.zeros(len(X_region)),
                            beta=pl.zeros(len(X_study)),
                            gamma=-5*pl.ones(len(est_mesh)),
                            sigma_alpha=[1.],
                            sigma_beta=[1.],
                            sigma_gamma=[10.],
                            delta=1000.,
                            sigma_delta=1.
                            )
    # cause-specific-mortality is a lower bound on p*f
    @mc.deterministic(name=key % 'pf')
    def pf(p=p, f=f):
        return (p+dismod3.settings.NEARLY_ZERO)*f

    data = [d for d in data_list if d['data_type'] == 'prevalence x excess-mortality data']
    lower_bound_data = [d for d in data_list if d['data_type'] == 'cause-specific mortality data']

    # make delta huge, sigma_delta tiny for pf (could hack sigma_delta = 0 to mean very precise)
    blank_prior_dict['delta'] = 1.e7
    blank_prior_dict['sigma_delta'] = 0.
    vars[key % 'prevalence_x_excess-mortality'] = neg_binom_model.setup(dm, key % 'pf', rate_stoch=pf, data_list=data, lower_bound_data=lower_bound_data, emp_prior=blank_prior_dict)
    blank_prior_dict['delta'] = 1000. # reset sigma_delta
    blank_prior_dict['sigma_delta'] = 1.    
    # m = m_all_cause - f * p
    @mc.deterministic(name=key % 'm_background')
    def m(SCpm=SCpm, param_mesh=dm.get_param_age_mesh(), est_mesh=dm.get_estimate_age_mesh()):
        return dismod3.utils.interpolate(param_mesh,  SCpm[3,:], est_mesh)
    vars[key % 'm_background'] = m

    # m_with = m + f
    @mc.deterministic(name=key % 'm_with')
    def m_with(m=m, f=f):
        return m + f
    data = [d for d in data_list if d['data_type'] == 'mortality data']
    # TODO: test this
    #prior_dict = dm.get_empirical_prior('excess-mortality')  # TODO:  make separate prior for with-condition mortality
    vars[key % 'mortality'] = neg_binom_model.setup(dm, key % 'm_with', data, m_with, emp_prior=blank_prior_dict)

    # mortality rate ratio = mortality with condition / mortality without
    @mc.deterministic(name=key % 'RR')
    def RR(m=m, m_with=m_with):
        return m_with / m
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
        pr_not_exit = pl.exp(-hazard)
        X = pl.empty(len(hazard))
        X[-1] = 1 / hazard[-1]
        for i in reversed(range(len(X)-1)):
            X[i] = pr_not_exit[i] * (X[i+1] + 1) + 1 / hazard[i] * (1 - pr_not_exit[i]) - pr_not_exit[i]
        return X
    data = [d for d in data_list if d['data_type'] == 'duration data']
    vars[key % 'duration'] = normal_model.setup(dm, key % 'duration', data, X)

    # YLD[a] = disability weight * i[a] * X[a] * regional_population[a]
    @mc.deterministic(name=key % 'i*X')
    def iX(i=i, X=X, p=p, pop=neg_binom_model.regional_population(key)):
        birth_yld = pl.zeros_like(p)
        birth_yld[0] = p[0] * pop[0]

        return i * X * (1-p) * pop + birth_yld
    
    vars[key % 'incidence_x_duration'] = {'rate_stoch': iX}

    return vars


