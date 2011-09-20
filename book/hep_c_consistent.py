import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}
models = {}

for ii in range(2):
    ### @export 'load model'
    dm = dismod3.load_disease_model(16391)

    ### @export 'set expert priors'
    dm.params['global_priors']['level_value']['incidence'] = dict(value=0., age_before=1., age_after=99)

    dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Slightly'
    dm.params['global_priors']['heterogeneity']['prevalence'] = 'Slightly'
    dm.params['global_priors']['level_value']['prevalence'] = dict(value=0., age_before=0, age_after=100)
    dm.params['global_priors']['level_bounds']['prevalence'] = dict(lower=0., upper =.05)
    dm.params['global_priors']['increasing']['prevalence'] = dict(age_start=0, age_end=0)
    dm.params['global_priors']['decreasing']['prevalence'] = dict(age_start=100, age_end=100)
    dm.params['sex_effect_prevalence'] = dict(mean=1, upper_ci=1.0001, lower_ci=.9999)
    dm.params['time_effect_prevalence'] = dict(mean=1, upper_ci=1.0001, lower_ci=.9999)
    dm.params['region_effect_prevalence'] = dict(std=.0001)
    dm.params['covariates']['Study_level']['bias']['rate']['value'] = 0
    for cv in dm.params['covariates']['Country_level']:
        dm.params['covariates']['Country_level'][cv]['rate']['value'] = 0

    # set bounds on remission and excess-mortality in the second time through
    if ii == 0:
        dm.params['global_priors']['level_bounds']['remission'] = dict(lower=0., upper =.05)
        dm.params['global_priors']['level_bounds']['excess_mortality'] = dict(lower=0., upper =.05)
        dm.params['global_priors']['level_bounds']['relative_risk'] = dict(lower=1., upper=1000.)
    if ii == 1:
        dm.params['global_priors']['level_bounds']['remission'] = dict(lower=0., upper =.01)
        dm.params['global_priors']['level_bounds']['excess_mortality'] = dict(lower=0., upper =.01)
        dm.params['global_priors']['level_bounds']['relative_risk'] = dict(lower=1., upper=1000.)

    ### @export 'initialize model data'
    region = 'north_america_high_income'
    year = 1990
    dm.data = [d for d in dm.data if dm.relevant_to(d, 'all', region, year, 'all')]

    # fit model
    dm.clear_fit()
    dm.clear_empirical_prior()
    dismod3.neg_binom_model.covariate_hash = {}

    import fit_world
    fit_world.fit_world(dm, generate_diagnostic_plots=False, store_results=False, map_only=False)
    models[ii] = dm
    #results[ii] = dict(rate_stoch=dm.vars['prevalence+world+all+all']['rate_stoch'].stats(), dispersion=dm.vars['prevalence+world+all+all']['dispersion'].stats())

    ### @export 'save'
    for d in dm.data:
        d['sex'] = 'male'  # otherwise tile plot shows it twice


    reload(book_graphics)
    book_graphics.plot_age_patterns(dm, region='world', year='all', sex='all',
                                    rate_types='remission incidence prevalence'.split())   
    pl.subplot(1,3,1)
    pl.yticks([0, .04, .08], [0,4,8])
    pl.ylabel('Rate (Per 100 PY)')

    #pl.subplot(1,4,2)
    #pl.yticks([0, .04, .08], [0,4,8])

    pl.subplot(1,3,2)
    pl.yticks([0, .0025, .0051], [0,.25,.5])

    pl.subplot(1,3,3)
    pl.yticks([0, .025, .05], [0,2.5,5])

    pl.savefig('hep_c-consistent%d.pdf' % ii)

pl.show()

book_graphics.save_json('hep_c.json', results)
