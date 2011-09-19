import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}
models = {}

sigma_vals = [.01, .1, 1.]

for sigma_sex in sigma_vals:
    for sigma_year in sigma_vals:
        ### @export 'load model'
        dm = dismod3.load_disease_model(16391)

        ### @export 'set expert priors'
        dm.set_param_age_mesh(pl.arange(0,101,10))
        dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Moderately'
        dm.params['global_priors']['heterogeneity']['prevalence'] = 'Slightly'

        dm.params['global_priors']['level_value']['prevalence'] = dict(value=0., age_before=0, age_after=100)
        dm.params['global_priors']['level_bounds']['prevalence'] = dict(lower=0., upper =.1)
        dm.params['global_priors']['increasing']['prevalence'] = dict(age_start=0, age_end=0)
        dm.params['global_priors']['decreasing']['prevalence'] = dict(age_start=100, age_end=100)
        dm.params['sex_effect_prevalence'] = dict(mean=1, upper_ci=pl.exp(sigma_sex*1.96), lower_ci=pl.exp(-sigma_sex*1.96))
        dm.params['time_effect_prevalence'] = dict(mean=1, upper_ci=pl.exp(sigma_year*1.96), lower_ci=pl.exp(-sigma_year*1.96))
        dm.params['covariates']['Study_level']['bias']['rate']['value'] = 1
        for cv in dm.params['covariates']['Country_level']:
            dm.params['covariates']['Country_level'][cv]['rate']['value'] = 0


        ### @export 'initialize model data'
        region = 'europe_eastern'
        dm.data = [d for d in dm.data if dm.relevant_to(d, 'prevalence', region, 'all', 'all')]

        # fit model
        dm.clear_fit()
        dm.clear_empirical_prior()
        dismod3.neg_binom_model.covariate_hash = {}

        dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence')
        models[sigma_sex, sigma_year] = dm

### @export 'save'
pl.figure(**book_graphics.half_page_params)

for row, sigma_sex in enumerate(sigma_vals):
    for col, sigma_year in enumerate(sigma_vals):
        index = row*3+col
        pl.subplot(3,3,1+index)

        if col == 0:
            if row == 1:
                pl.ylabel('Prevalence (Per 100)')
            pl.yticks([0, .02, .04, .06], [0, 2, 4, 6]) 
        else:
            pl.yticks([0, .02, .04, .06], ['', '', '', '']) 

        if row == 2:
            pl.xlabel('Age (Years)')
            pl.xticks([0,25,50,75])
        else:
            pl.xticks([0,25,50,75], ['','','',''])

        dismod3.plotting.plot_intervals(dm, [d for d in dm.data if dm.relevant_to(d, 'prevalence', region, 'all', 'all')],
                                        color='black', print_sample_size=False, alpha=1., plot_error_bars=False,
                                        linewidth=2)
        for r in models[sigma_sex, sigma_year].vars['rate_stoch'].trace()[::10]:
            pl.step(range(101), r, '-', color='grey', linewidth=2, zorder=-100)
        pl.step(range(101), models[sigma_sex, sigma_year].vars['rate_stoch'].stats()['quantiles'][50],
                linewidth=3, color='white')
        pl.step(range(101), models[sigma_sex, sigma_year].vars['rate_stoch'].stats()['quantiles'][50],
                linewidth=1, color='black')

        pl.axis([0, 100, 0, .075])
        pl.text(5, .07, '(%s)'%('abcdefghi'[index]), va='top', ha='left')

pl.subplots_adjust(hspace=0, wspace=0, bottom=.1, left=.06, right=.99, top=.97)
pl.savefig('hep_c-year_sex.pdf')
pl.savefig('hep_c-year_sex.png')
