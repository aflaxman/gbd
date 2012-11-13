import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

region = 'asia_southeast'
def initialize_model():
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
    dm.params['covariates']['Study_level']['bias']['rate']['value'] = 1
    for cv in dm.params['covariates']['Country_level']:
        dm.params['covariates']['Country_level'][cv]['rate']['value'] = 0


    ### @export 'initialize model data'
    dm.data = [d for d in dm.data if dm.relevant_to(d, 'prevalence', region, 'all', 'all')]

    # fit model
    dm.clear_fit()
    dm.clear_empirical_prior()
    dismod3.neg_binom_model.covariate_hash = {}
    return dm

models = []

models.append([])
for level in ['Slightly', 'Moderately', 'Very']:
    dm = initialize_model()
    dm.params['global_priors']['smoothness']['prevalence']['amount'] = level
    dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence', map_only=True)
    models[-1].append(dm)

models.append([])
for level in [.01, .1, 1.]:
    dm = initialize_model()
    dm.params['sex_effect_prevalence'] = dict(mean=1, upper_ci=pl.exp(level *1.96), lower_ci=pl.exp(-level*1.96))
    dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence', map_only=True)
    models[-1].append(dm)
    
models.append([])
for level in [.01, .1, 1.]:
    dm = initialize_model()
    dm.params['time_effect_prevalence'] = dict(mean=1, upper_ci=pl.exp(level *1.96), lower_ci=pl.exp(-level*1.96))
    dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence', map_only=True)
    models[-1].append(dm)



models.append([])
for level in [.01, .1, 1.]:
    dm = initialize_model()
    dm.params['region_effect_prevalence'] = dict(std=level)
    dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence', map_only=True)
    models[-1].append(dm)

models.append([])
for level in [.1, 1, 10.]:
    dm = initialize_model()
    dm.params['beta_effect_prevalence'] = dict(mean=[0.], std=[level])
    dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence', map_only=True)
    models[-1].append(dm)

models.append([])
for level in [.2, 2., 20.]:
    dm = initialize_model()
    dm.params['gamma_effect_prevalence'] = dict(mean=list(pl.zeros_like(dm.get_estimate_age_mesh())),
                                                std=list(level*pl.ones_like(dm.get_estimate_age_mesh())))
    dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence', map_only=True)
    models[-1].append(dm)



# this should change uncertainty, although it turns out not to change levels
models.append([])
for level in [.025, .25, 2.5]:
    dm = initialize_model()
    dm.params['delta_effect_prevalence'] = dict(mean=3.,
                                                std=level)
    dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence', map_only=True)
    models[-1].append(dm)

models.append([])
for level in [1, 2, 3]:
    dm = initialize_model()
    dm.params['delta_effect_prevalence'] = dict(mean=level,
                                                std=.25)
    dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence', map_only=True)
    models[-1].append(dm)

models.append([])
for level in [.025, .25, 2.5]:
    dm = initialize_model()
    dm.params['zeta_effect_prevalence'] = dict(mean=0.,
                                                std=level)
    dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence', map_only=True)
    models[-1].append(dm)


### @export 'save'
pl.figure(**book_graphics.half_page_params)

color=['red', 'green', 'blue']
for row in range(3):
    for col in range(3):
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
                                        linewidth=2, zorder=-10)
        for level in range(3):
            dm = models[index][level]
            pl.step(range(101), dm.vars['rate_stoch'].value,
                    linewidth=3, color='white', zorder=-1)
            pl.step(range(101), dm.vars['rate_stoch'].value,
                    linewidth=1, color=color[level], alpha=.5)

        pl.axis([0, 100, 0, .075])
        pl.text(5, .07, '(%s)'%('abcdefghi'[index]), va='top', ha='left')

pl.subplots_adjust(hspace=0, wspace=0, bottom=.1, left=.06, right=.99, top=.97)
pl.savefig('hep_c-sensitivity.pdf')
pl.savefig('hep_c-sensitivity.png')
