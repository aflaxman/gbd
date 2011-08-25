import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}
models = {}

smoothness = ['No Prior', 'Slightly', 'Moderately', 'Very']
linestyle = dict(zip(smoothness, ['steps-mid-', 'steps-mid:', 'steps-mid--', 'steps-mid-.']))

for smooth_i in smoothness:
    ### @export 'load model'
    dm = dismod3.load_disease_model(16391)

    ### @export 'set expert priors'
    dm.params['global_priors']['smoothness']['prevalence']['amount'] = smooth_i
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


    ### @export 'initialize model data'
    region = 'north_america_high_income'
    year = 1990
    dm.data = [d for d in dm.data if dm.relevant_to(d, 'prevalence', region, year, 'all')]

    # fit model
    dm.clear_fit()
    dm.clear_empirical_prior()
    dismod3.neg_binom_model.covariate_hash = {}

    dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence')
    models[smooth_i] = dm
    results[smooth_i] = dict(rate_stoch=dm.vars['rate_stoch'].stats(), dispersion=dm.vars['dispersion'].stats())

### @export 'save'
for d in dm.data:
    d['sex'] = 'male'  # otherwise tile plot shows it twice

dismod3.plotting.tile_plot_disease_model(dm, dismod3.utils.gbd_keys(['prevalence'], [region], [str(year)], ['all']),
                                         plot_prior_flag=False, print_sample_size=False, plot_error_bars=False)

for ii, smooth_i in enumerate(smoothness):
    pl.plot(pl.arange(101)+ii, models[smooth_i].vars['rate_stoch'].stats()['mean'],
            linewidth=3, color='white', linestyle='steps-mid')
    pl.plot(pl.arange(101)+ii, models[smooth_i].vars['rate_stoch'].stats()['mean'],
            linewidth=1, color='black', linestyle=linestyle[smooth_i],
            label=smooth_i)

pl.legend(loc='upper left', title='Smoothing Prior', fancybox=True, shadow=True)

pl.axis([0, 100, 0, .075])
pl.title('')
pl.ylabel('Prevalence (Per 1)')
pl.xlabel('Age (Years)')
pl.savefig('hep_c-smoothing.pdf')
