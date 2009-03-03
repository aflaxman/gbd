import inspect

import numpy as np
import pymc as mc

import probabilistic_utils
import single_binomial_rate as rate_model

MCMC_PARAMS = {
    'most accurate': [500, 20, 10000],
    'testing fast': [5, 5, 5],
    }

MAP_PARAMS = {
    'most accurate': [ 500, 'fmin_powell'],
    'try powells method': [ 100, 'fmin_powell'],
    'testing fast': [ 1, 'fmin' ],
    }

def map_fit(dm, speed='most accurate'):
    """
    tuck the pymc vars into the disease_model.vars, if they don't
    exist and then fit them with map
    """
    setup_disease_model(dm)
    map = mc.MAP(dm.vars)
    
    print "searching for maximum likelihood point estimate (%s)" % speed
    iterlim, method = MAP_PARAMS[speed]
    map.fit(verbose=10, iterlim=iterlim, method=method)

    probabilistic_utils.save_map(dm.i)
    probabilistic_utils.save_map(dm.r)
    probabilistic_utils.save_map(dm.f)
    probabilistic_utils.save_map(dm.p)

def mcmc_fit(dm, speed='most accurate'):
    #speed = 'testing fast'
    map_fit(dm, speed)

    # clear any data from previous mcmc fit
    # TODO: refactor this, or really switch to an open notebook
    # approach, where we save _all_ fits
    for rf in [dm.i, dm.r, dm.p, dm.f]:
        rf.fit.pop('mcmc_mean', '')
        rf.fit.pop('mcmc_median', '')
        rf.fit.pop('mcmc_upper_cl', '')
        rf.fit.pop('mcmc_lower_cl', '')
        rf.save()
    
    trace_len, thin, burn = MCMC_PARAMS[speed]
    mcmc = mc.MCMC(dm.vars)

    # TODO: refactor the step method setup code into the rate_model
    for rf in [dm.i, dm.r, dm.p, dm.f]:
        try:
            mcmc.use_step_method(mc.AdaptiveMetropolis, rf.vars['beta_binom_stochs'], verbose=0)
        except (KeyError, ValueError):
            pass

    mcmc.sample(trace_len*thin+burn, burn, thin, verbose=1)

    probabilistic_utils.save_mcmc(dm.i)
    probabilistic_utils.save_mcmc(dm.r)
    probabilistic_utils.save_mcmc(dm.f)
    probabilistic_utils.save_mcmc(dm.p)

def initialized_rate_vars(rf, rate_stoch=None):
    # store the rate model code in the asrf for future reference
    try:
        rf.fit['rate_model'] = inspect.getsource(rate_model)
    except:
        rf.fit['rate_model'] = 'WARNING: could not find source code'

    rf.fit['out_age_mesh'] = range(probabilistic_utils.MAX_AGE)

    # do normal approximation first, to generate a good starting point
    M,C = probabilistic_utils.normal_approx(rf)
    rate_model.setup_rate_model(rf, rate_stoch)
    return probabilistic_utils.flatten(rf.vars.values()), rf.rate_stoch
    

#################### Code for generating a hierarchical bayesian model of the asrf interactions
def setup_disease_model(dm):
    """
    Generate a list of the PyMC variables for the generic disease
    model, and store it as dm.vars

    See comments in the code for exact details on the model.
    """
    dm.vars = {}
    out_age_mesh = range(probabilistic_utils.MAX_AGE)

    dm.vars['i'], i = initialized_rate_vars(dm.i_in())
    dm.vars['r'], r = initialized_rate_vars(dm.r_in())
    dm.vars['f'], f = initialized_rate_vars(dm.f_in())
    # TODO: create m, from all-cause mortality
    m = np.zeros(len(out_age_mesh))
    
    # TODO: make error in C_0 a semi-informative stochastic variable
    logit_C_0 = mc.Normal('logit(C_0)', 0., 1.e-2)
    @mc.deterministic
    def C_0(logit_C_0=logit_C_0):
        return mc.invlogit(logit_C_0)
    
    @mc.deterministic
    def S_0(C_0=C_0):
        return max(0.0, 1.0 - C_0)
    dm.vars['bins'] = [S_0, C_0, logit_C_0]
    
    # iterative solution to difference equations to obtain bin sizes for all ages
    @mc.deterministic
    def S_C_D_M(S_0=S_0, C_0=C_0, i=i, r=r, f=f, m=m):
        S = np.zeros(len(out_age_mesh))
        C = np.zeros(len(out_age_mesh))
        D = np.zeros(len(out_age_mesh))
        M = np.zeros(len(out_age_mesh))
        
        S[0] = S_0
        C[0] = C_0
        D[0] = 0.0
        M[0] = 0.0
        
        for a in range(len(out_age_mesh)-1):
            S[a+1] = S[a]*(1-i[a]-m[a]) + C[a]*r[a]
            C[a+1] = S[a]*i[a]          + C[a]*(1-r[a]-m[a]-f[a])
            D[a+1] =                      C[a]*f[a]               + D[a]
            M[a+1] = S[a]*m[a]          + C[a]*m[a]                      + M[a]
                
        return S,C,D,M
    dm.vars['bins'] += [S_C_D_M]

    # prevalence = # with condition / (# with condition + # without)
    @mc.deterministic
    def p(S_C_D_M=S_C_D_M, tau_p=1./.01**2):
        S,C,D,M = S_C_D_M
        return probabilistic_utils.trim(C/(S+C),
                                        probabilistic_utils.NEARLY_ZERO,
                                        1. - probabilistic_utils.NEARLY_ZERO)
    dm.vars['p'], p = initialized_rate_vars(dm.p_in(), rate_stoch=p)

