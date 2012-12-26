import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

region = 'latin_america_southern'
year = 'all'
def initialize_model():
    ### @export 'load model'
    dm = dismod3.load_disease_model(16391)

    ### @export 'set expert priors'
    dm.set_param_age_mesh(pl.arange(0,101,10))
    dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Slightly'
    dm.params['global_priors']['heterogeneity']['prevalence'] = 'Slightly'

    dm.params['global_priors']['level_value']['prevalence'] = dict(value=0., age_before=0, age_after=100)
    dm.params['global_priors']['level_bounds']['prevalence'] = dict(lower=0., upper =.1)
    dm.params['global_priors']['increasing']['prevalence'] = dict(age_start=0, age_end=0)
    dm.params['global_priors']['decreasing']['prevalence'] = dict(age_start=100, age_end=100)
    dm.params['covariates']['Study_level']['bias']['rate']['value'] = 1
    for cv in dm.params['covariates']['Country_level']:
        dm.params['covariates']['Country_level'][cv]['rate']['value'] = 0


    ### @export 'initialize model data'
    dm.data = [d for d in dm.data if dm.relevant_to(d, 'prevalence', region, year, 'all')]

    # fit model
    dm.clear_fit()
    dm.clear_empirical_prior()
    dismod3.neg_binom_model.covariate_hash = {}
    return dm

models = []

for mean in [1., 2., 3.]:
    for std in [.025, .25, 2.5]:
        dm = initialize_model()
        dm.params['delta_effect_prevalence'] = dict(mean=mean,
                                                    std=std)
        dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence', store_results=False)
        models.append(dm)

### @export 'save'
### @export 'save'
pl.figure(**book_graphics.half_page_params)
for ii, dm in enumerate(models):
    pl.subplot(3, 3, ii+1)
    dismod3.plotting.plot_intervals(dm, [d for d in dm.data if dm.relevant_to(d, 'prevalence', region, year, 'all')],
                                    color='black', print_sample_size=False, alpha=1., plot_error_bars=False,
                                    linewidth=2)
    for r in dm.vars['rate_stoch'].trace()[::10]:
        pl.step(range(101), r, '-', color='grey', linewidth=2, zorder=-100)
    pl.step(range(101), dm.vars['rate_stoch'].stats()['quantiles'][50],
            linewidth=3, color='white')
    pl.step(range(101), dm.vars['rate_stoch'].stats()['quantiles'][50],
            linewidth=1, color='black')
    pl.axis([0, 100, 0, .075])
    pl.title('')
    if (ii/3) < 2:
        pl.xticks([0,25,50,75], ['', '', '', ''])
    else:
        pl.xticks([0,25,50,75])
        pl.xlabel('Age (Years)')

    if (ii%3) > 0:
        pl.yticks([0, .02, .04, .06], ['', '', '', '']) 
    else:
        pl.yticks([0, .02, .04, .06], [0, 2, 4, 6]) 
        if (ii/3) == 1:
            pl.ylabel('Prevalence (Per 100)')
    pl.text(5, .07, '(%s)'% ('abcdefghi'[ii]), va='top', ha='left')

pl.subplots_adjust(hspace=0, wspace=0, bottom=.1, left=.06, right=.99, top=.97)
pl.savefig('hep_c-smoothing.pdf')
