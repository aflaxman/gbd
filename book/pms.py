import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}

### @export 'data'
# load model
dm = dismod3.load_disease_model(16370)

# set expert priors and other model parameters
dm.set_param_age_mesh([0, 15, 20, 25, 35, 45, 50, 55, 100])
#dm.set_param_age_mesh([0, 15, 20, 25, 30, 35, 40, 45, 50, 55, 100])
#dm.set_param_age_mesh([0, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 100])

dm.params['global_priors']['level_value']['incidence']['age_before'] = 15
dm.params['global_priors']['level_value']['incidence']['age_after'] = 50
#dm.params['global_priors']['smoothness']['incidence']['age_start'] = 15

dm.params['global_priors']['level_value']['remission']['age_before'] = 40
dm.params['global_priors']['level_bounds']['remission']['upper'] = 10.

dm.params['global_priors']['level_value']['excess_mortality']['age_before'] = 101

dm.params['global_priors']['level_value']['prevalence']['age_before'] = 15
dm.params['global_priors']['level_value']['prevalence']['age_after'] = 50

dm.params['covariates']['Country_level']['LDI_id_Updated_7July2011']['rate']['value'] = 0


# clear any fit and priors
dm.clear_fit()
dm.clear_empirical_prior()
dismod3.neg_binom_model.covariate_hash = {}

# initialize model data
prev_data = [d for d in dm.data if d['data_type'] == 'prevalence data']
r = pl.array([dm.value_per_1(s) for s in prev_data])
min_rate_per_100 = '%d' % round(r.min()*100)
max_rate_per_100 = '%d' % round(r.max()*100)

import fit_world
#fit_world.fit_world(dm)

import fit_posterior
#fit_posterior.fit_posterior(dm, 'asia_south', 'female', '2005')
#fit_posterior.fit_posterior(dm, 'europe_western', 'female', '2005')
fit_posterior.fit_posterior(dm, 'north_america_high_income', 'female', '2005')

pl.figure(**book_graphics.quarter_page_params)
pl.subplot(1,2,1)
dismod3.plotting.plot_intervals(dm, [d for d in dm.data if dm.relevant_to(d, 'prevalence', 'all', 'all', 'all')],
                                color='black', print_sample_size=False, alpha=.75, plot_error_bars=False,
                                linewidth=1)
pl.axis([10,60,-.01,1])
pl.yticks([0,.25,.5,.75])
pl.ylabel('Prevalence (per 1)')
pl.xlabel('Age (years)')
pl.title('a) All data')


pl.subplot(1,2,2)
# FIXME: fit_posterior removed all data
dismod3.plotting.plot_intervals(dm, [d for d in dm.data if dm.relevant_to(d, 'prevalence', 'north_america_high_income', 'all', 'all')],
                                color='black', print_sample_size=False, alpha=.75, plot_error_bars=False,
                                linewidth=1)
for r in dm.vars['prevalence+north_america_high_income+2005+female']['rate_stoch'].trace():
    plot(dm.get_estimate_age_mesh(), r, '-', color='grey', linewidth=2, zorder=-100)

r = dm.vars['prevalence+north_america_high_income+2005+female']['rate_stoch'].stats()['quantiles'][50]
plot(dm.get_estimate_age_mesh(), r, 'w-', linewidth=2, zorder=99)
plot(dm.get_estimate_age_mesh(), r, 'k-', linewidth=1, zorder=100)

pl.axis([10,60,-.01,1])
pl.yticks([0,.25,.5,.75])
pl.ylabel('Prevalence (per 1)')
pl.xlabel('Age (years)')
pl.title('b) North America, High Income')

pl.subplots_adjust(bottom=.15, wspace=.25, top=.85, left=.1, right=.95)

pl.savefig('pms-prev.pdf')
pl.show()

### @export 'save'
book_graphics.save_json('pms.json', vars())
