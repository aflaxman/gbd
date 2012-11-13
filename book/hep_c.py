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
    dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Slightly'
    dm.params['global_priors']['heterogeneity']['prevalence'] = 'Slightly'
    dm.params['global_priors']['level_value']['prevalence'] = dict(value=0., age_before=0, age_after=100)
    dm.params['global_priors']['level_bounds']['prevalence'] = dict(lower=0., upper =.1)
    dm.params['global_priors']['increasing']['prevalence'] = dict(age_start=0, age_end=10)
    dm.params['global_priors']['decreasing']['prevalence'] = dict(age_start=80, age_end=100)
    dm.params['sex_effect_prevalence'] = dict(mean=1, upper_ci=1.0001, lower_ci=.9999)
    dm.params['time_effect_prevalence'] = dict(mean=1, upper_ci=1.0001, lower_ci=.9999)
    dm.params['covariates']['Study_level']['bias']['rate']['value'] = 0
    for cv in dm.params['covariates']['Country_level']:
        dm.params['covariates']['Country_level'][cv]['rate']['value'] = 0

    # TODO: set bounds on remission and excess-mortality in the second time through

    ### @export 'initialize model data'
    region = 'north_america_high_income'
    year = 1990
    dm.data = [d for d in dm.data if dm.relevant_to(d, 'prevalence', region, year, 'all')]

    # fit model
    dm.clear_fit()
    dm.clear_empirical_prior()
    dismod3.neg_binom_model.covariate_hash = {}

    import fit_world
    fit_world.fit_world(dm)
    models[ii] = dm
    results[ii] = dict(rate_stoch=dm.vars['prevalence+world+all+all']['rate_stoch'].stats(), dispersion=dm.vars['prevalence+world+all+all']['dispersion'].stats())

    ### @export 'save'
    for d in dm.data:
        d['sex'] = 'male'  # otherwise tile plot shows it twice

    dismod3.plotting.tile_plot_disease_model(dm, dismod3.utils.gbd_keys(['incidence', 'remission', 'excess-mortality', 'prevalence'], ['world'], ['all'], ['all']),
                                             plot_prior=False, print_sample_size=False)
pl.show()

book_graphics.save_json('hep_c.json', vars())
