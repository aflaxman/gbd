import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}
models = {}

smoothness = ['Slightly', 'Moderately', 'Very']

for smooth_i in smoothness:
    ### @export 'load model'
    dm = dismod3.load_disease_model(16391)

    ### @export 'set expert priors'
    dm.params['global_priors']['smoothness']['prevalence']['amount'] = smooth_i
    dm.params['global_priors']['heterogeneity']['prevalence'] = 'Slightly'
    dm.params['global_priors']['level_value']['prevalence'] = dict(value=0., age_before=0, age_after=100)
    dm.params['global_priors']['level_bounds']['prevalence'] = dict(lower=0., upper =.1)
    dm.params['global_priors']['increasing']['prevalence'] = dict(age_start=0, age_end=0)
    dm.params['global_priors']['decreasing']['prevalence'] = dict(age_start=100, age_end=100)
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
pl.figure(**book_graphics.quarter_page_params)
pl.subplot(1,4,1)
ax=dismod3.plotting.plot_intervals(dm, [d for d in dm.data if dm.relevant_to(d, 'prevalence', region, year, 'all')],
                                color='black', print_sample_size=False, alpha=1., plot_error_bars=True,
                                linewidth=2)
pl.ylabel('Prevalence (Per 100)')
pl.title('')
pl.xlabel('Age (Years)')
pl.axis([0, 100, 0, .075])
pl.xticks([0,25,50,75])
pl.yticks([0, .02, .04, .06], [0, 2, 4, 6]) 
pl.text(5, .07, 'a)', va='top', ha='left')

for ii, smooth_i in enumerate(smoothness):
    pl.subplot(1,4,ii+2)
    dismod3.plotting.plot_intervals(dm, [d for d in dm.data if dm.relevant_to(d, 'prevalence', region, year, 'all')],
                                    color='black', print_sample_size=False, alpha=1., plot_error_bars=False,
                                    linewidth=2)
    for r in models[smooth_i].vars['rate_stoch'].trace():
        pl.step(range(101), r, '-', color='grey', linewidth=2, zorder=-100)
    pl.step(range(101), models[smooth_i].vars['rate_stoch'].stats()['quantiles'][50],
            linewidth=3, color='white')
    pl.step(range(101), models[smooth_i].vars['rate_stoch'].stats()['quantiles'][50],
            linewidth=1, color='black')
    pl.axis([0, 100, 0, .075])
    pl.title('')
    pl.xlabel('Age (Years)')
    pl.xticks([0,25,50,75])
    pl.yticks([0, .02, .04, .06], ['', '', '', '']) 
    pl.text(5, .07, '%s)'% ('bcd'[ii]), va='top', ha='left')

pl.subplots_adjust(wspace=0, bottom=.15, left=.06, right=.99, top=.97)
pl.savefig('hep_c-smoothing.pdf')
