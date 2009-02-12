import inspect

import numpy as np
import pymc as mc

import probabilistic_utils
import rate_single_binomial as rate_model

def map_fit(asrfs):
    all_vars = []
    rate_vars = {}
    rate_stochs = []

    for asrf in asrfs:
        # store the rate model code in the asrf for future reference
        asrf.fit['rate_model'] = inspect.getsource(rate_model)
        asrf.fit['out_age_mesh'] = range(probabilistic_utils.MAX_AGE)

        # do normal approximation first, to generate a good starting point
        M,C = probabilistic_utils.normal_approx(asrf)

        # define the model variables
        vars = rate_model.model_vars(asrf)
        rate_vars[asrf.id] = vars
        rate_stochs += [ vars['asrf_%d'%asrf.id] ]
        all_vars += vars.values()

    @mc.potential
    def smooth_across_regions(rate_list=rate_stochs):
        logp = 0.
        for ii in range(len(rate_list)):
            for jj in range(ii+1, len(rate_list)):
                logp += mc.normal_like(np.diff(np.log(rate_list[ii]))-np.diff(np.log(rate_list[jj])), 0., 1./(.1)**2)
        return logp
    all_vars += [ smooth_across_regions ]

    map = mc.MAP(all_vars)

    print "searching for maximum likelihood point estimate"
    iterlim = 500
    method = 'fmin_powell'
    map.fit(verbose=10, iterlim=iterlim, method=method)

    for asrf in asrfs:
        rate_model.save_map(rate_vars[asrf.id], asrf)

    return all_vars, rate_vars

def mcmc_fit(asrfs):
    all_vars, rate_vars = map_fit(asrfs)
    
    trace_len, thin, burn = 500, 10, 5000
    mcmc = mc.MCMC(all_vars)
    mcmc.sample(trace_len*thin+burn, burn, thin, verbose=1)

    for asrf in asrfs:
        rate_model.save_mcmc(rate_vars[asrf.id], asrf)

