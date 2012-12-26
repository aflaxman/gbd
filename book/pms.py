import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

faster_run_flag = False

results = {}

### @export 'data'
# load model
dm = dismod3.load_disease_model(16370)

# set expert priors and other model parameters
dm.set_param_age_mesh([0, 15, 17, 20, 25, 35, 45, 50, 55, 100])

dm.params['global_priors']['level_value']['incidence']['age_before'] = 15
dm.params['global_priors']['level_value']['incidence']['age_after'] = 50
dm.params['global_priors']['smoothness']['incidence']['age_start'] = 15

dm.params['global_priors']['level_value']['remission']['age_before'] = 40
dm.params['global_priors']['level_bounds']['remission']['upper'] = 10.

dm.params['global_priors']['level_value']['excess_mortality']['age_before'] = 101

dm.params['global_priors']['level_value']['prevalence']['age_before'] = 15
dm.params['global_priors']['level_value']['prevalence']['age_after'] = 50
dm.params['global_priors']['smoothness']['prevalence']['age_start'] = 15

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
median_rate_per_100 = '%d' % round(pl.median(r*100))
regions = pl.array([d['gbd_region'] for d in prev_data])
num_regions = len(pl.unique(regions))

import fit_world
#fit_world.fit_world(dm)
#dm.data = prev_data # put data back in

import fit_posterior
region = 'north_america_high_income'
sex = 'female'
year='2005'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=faster_run_flag, store_results=False)
dm.data = prev_data # put data back in

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
dismod3.plotting.plot_intervals(dm, [d for d in dm.data if dm.relevant_to(d, 'prevalence', region, year, sex)],
                                color='black', print_sample_size=False, alpha=.75, plot_error_bars=False,
                                linewidth=1)
book_graphics.plot_rate(dm, dismod3.utils.gbd_key_for('prevalence', region, year, sex), linestyle='-')

pl.axis([10,60,-.01,1])
pl.yticks([0,.25,.5,.75])
pl.ylabel('Prevalence (per 1)')
pl.xlabel('Age (years)')
pl.title('b) North America, High Income, %s'%year)

pl.subplots_adjust(bottom=.15, wspace=.25, top=.85, left=.1, right=.95)

pl.savefig('pms-prev.pdf')

book_graphics.plot_age_patterns(dm, region=region, year=year, sex=sex,
                                xticks=[15, 25, 35, 45, 55], rate_types='incidence remission prevalence'.split(),
                                yticks=dict(incidence=[0, .05, .1, .15, .2], remission=[0, .5, 1, 1.5, 2.], prevalence=[0, .25, .5, .75, 1]))

pl.savefig('pms-consistent.pdf')

### @export 'save'
book_graphics.save_json('pms.json', vars())
