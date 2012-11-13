import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import fit_posterior

import book_graphics
reload(book_graphics)

region = 'europe_western'
year = '2005'
sex = 'male'
ymax = .002

def initialize_model():
    ### @export 'load model'
    dm = dismod3.load_disease_model(19807)
    ### @export 'initialize model data'
    dm.params['global_priors']['level_bounds']['excess_mortality'] = dict(lower=.1, upper=100.)
    dm.params['global_priors']['increasing']['excess_mortality'] = dict(age_start=0, age_end=0)
    dm.params['global_priors']['level_bounds']['relative_risk'] = dict(lower=0., upper=10000.)

    for cv in dm.params['covariates']['Study_level']:
        dm.params['covariates']['Study_level'][cv]['rate']['value'] = 0
    for cv in dm.params['covariates']['Country_level']:
        dm.params['covariates']['Country_level'][cv]['rate']['value'] = 0

    level = .001
    dm.params['sex_effect_prevalence'] = dict(mean=1, upper_ci=pl.exp(level *1.96), lower_ci=pl.exp(-level*1.96))
    dm.params['time_effect_prevalence'] = dict(mean=1, upper_ci=pl.exp(level *1.96), lower_ci=pl.exp(-level*1.96))
    dm.params['region_effect_prevalence'] = dict(std=level)

    dm.clear_fit()
    dm.clear_empirical_prior()
    dismod3.neg_binom_model.covariate_hash = {}
    return dm


models = []

dm = initialize_model()
dm.description = 'As loaded, but only GBR data, no covariates'
dm.data = [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)] + \
    [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)] + \
    [d for d in dm.data if d.get('country_iso3_code') == 'GBR'] + \
    [d for d in dm.data if dm.relevant_to(d, 'excess-mortality', 'all', 'all', 'all')]
dm.params['global_priors']['increasing']['excess_mortality'] = dict(age_start=25, age_end=100)
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)

dm = initialize_model()
dm.description = 'Without increasing prior on excess-mortality'
dm.data = [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)] + \
    [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)] + \
    [d for d in dm.data if d.get('country_iso3_code') == 'GBR'] + \
    [d for d in dm.data if dm.relevant_to(d, 'excess-mortality', 'all', 'all', 'all')]
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)

dm = initialize_model()
dm.description = 'With lower bound of .2 on excess-mortality (to encourage good convergence)'
dm.data = [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)] + \
    [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)] + \
    [d for d in dm.data if d.get('country_iso3_code') == 'GBR'] + \
    [d for d in dm.data if dm.relevant_to(d, 'excess-mortality', 'all', 'all', 'all')]
dm.params['global_priors']['level_bounds']['excess_mortality'] = dict(lower=.2, upper=10.)
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)

dm = initialize_model()
dm.description = 'With excess-mortality data removed'
dm.data = [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)] + \
    [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)] + \
    [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
dm.params['global_priors']['level_bounds']['excess_mortality'] = dict(lower=.2, upper=10.)
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)

dm = initialize_model()
dm.description = 'With excess-mortality data and priors removed'
dm.data = [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)] + \
    [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)] + \
    [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
dm.params['global_priors']['level_bounds']['excess_mortality'] = dict(lower=0., upper=10.)
dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'No Prior'
dm.params['gamma_effect_excess-mortality'] = dict(mean=list(pl.zeros_like(dm.get_estimate_age_mesh())),
                                                  std=list(10.*pl.ones_like(dm.get_estimate_age_mesh())))
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)


dm = initialize_model()
dm.description = 'Without increasing prior on excess-mortality, but with "Very" smoothing'
dm.data = [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)] + \
    [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)] + \
    [d for d in dm.data if d.get('country_iso3_code') == 'GBR'] + \
    [d for d in dm.data if dm.relevant_to(d, 'excess-mortality', 'all', 'all', 'all')]
dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Very'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)

dm = initialize_model()
dm.description = 'Without increasing prior on excess-mortality, with "Very" smoothing, no excess-mortality data'
dm.data = [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)] + \
    [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)] + \
    [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Very'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)

dm = initialize_model()
dm.description = 'Without increasing prior on excess-mortality, with "No Prior" smoothing on incidence, no excess-mortality data'
dm.data = [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)] + \
    [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)] + \
    [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Slightly'
dm.params['global_priors']['smoothness']['incidence']['amount'] = 'No Prior'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)

dm = initialize_model()
dm.description = 'only prevalence data'
dm.data = [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)] + \
    [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Slightly'
dm.params['global_priors']['smoothness']['incidence']['amount'] = 'No Prior'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=False, store_results=False)
models.append(dm)

dm = initialize_model()
dm.description = 'up the prevalence sample size, to get fit to follow prev data when pf data is back in'
prev_data = [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
for d in prev_data:
    d['effective_sample_size'] = 1.e7
dm.params['delta_effect_prevalence'] = dict(mean=5.,
                                            std=.025)
dm.data = prev_data +[d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)] + \
    [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)]
dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Slightly'
dm.params['global_priors']['smoothness']['incidence']['amount'] = 'No Prior'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)

dm = initialize_model()
dm.description = 'shrink age intervals to make things simpler'
prev_data = [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
for d in prev_data:
    d['effective_sample_size'] = 1.e7
    d['age_end'] = d['age_start']
    d['age_weights'] = pl.array([1.])
dm.params['delta_effect_prevalence'] = dict(mean=5.,
                                            std=.025)
dm.data = prev_data +[d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)]
dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Moderately'
dm.params['global_priors']['smoothness']['incidence']['amount'] = 'No Prior'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)



dm = initialize_model()
dm.description = 'shrink age intervals to make things simpler, include pf and p'
prev_data = [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
for d in prev_data:
    d['effective_sample_size'] = 1.e7
dm.params['delta_effect_prevalence'] = dict(mean=5.,
                                            std=.025)

pf_data = [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)]
for d in pf_data + prev_data:
    d['age_end'] = d['age_start']
    d['age_weights'] = pl.array([1.])
dm.data = prev_data + pf_data + [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)]

dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Moderately'
dm.params['global_priors']['smoothness']['incidence']['amount'] = 'No Prior'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)



dm = initialize_model()
dm.description = 'slightly larger age intervals'
prev_data = [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
for d in prev_data:
    d['effective_sample_size'] = 1.e7
dm.params['delta_effect_prevalence'] = dict(mean=5.,
                                            std=.025)

pf_data = [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)]
for d in pf_data + prev_data:
    d['age_end'] = d['age_start']+1
    d['age_weights'] = pl.array([.5, .5])
dm.data = prev_data + pf_data + [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)]

dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Moderately'
dm.params['global_priors']['smoothness']['incidence']['amount'] = 'No Prior'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)


dm = initialize_model()
dm.description = 'return to good old age intervals'
prev_data = [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
for d in prev_data:
    d['effective_sample_size'] = 1.e7
dm.params['delta_effect_prevalence'] = dict(mean=5.,
                                            std=.025)

pf_data = [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)]
for d in pf_data + prev_data:
    d['effective_sample_size'] = 1.e7

dm.data = prev_data + pf_data + [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)]

dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Moderately'
dm.params['global_priors']['smoothness']['incidence']['amount'] = 'No Prior'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)


dm = initialize_model()
dm.description = 'return to old sample sizes, show that this is the problem'
prev_data = [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
dm.params['delta_effect_prevalence'] = dict(mean=5.,
                                            std=.025)
pf_data = [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)]
dm.data = prev_data + pf_data + [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)]

dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Moderately'
dm.params['global_priors']['smoothness']['incidence']['amount'] = 'No Prior'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)


dm = initialize_model()
dm.description = 'adjust sample sizes to explore more'
prev_data = [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
for d in prev_data:
    d['effective_sample_size'] = 1.e6
dm.params['delta_effect_prevalence'] = dict(mean=5.,
                                            std=.025)

pf_data = [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)]
for d in pf_data + prev_data:
    d['effective_sample_size'] = 1.e6

dm.data = prev_data + pf_data + [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)]

dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Moderately'
dm.params['global_priors']['smoothness']['incidence']['amount'] = 'No Prior'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)


dm = initialize_model()
dm.description = 'adjust sample sizes to explore more'
prev_data = [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
for d in prev_data:
    d['effective_sample_size'] = 1.e6
dm.params['delta_effect_prevalence'] = dict(mean=5.,
                                            std=.025)
dm.params['global_priors']['decreasing']['prevalence'] = dict(age_start=70, age_end=100)

pf_data = [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)]
for d in pf_data + prev_data:
    d['effective_sample_size'] = 1.e6

dm.data = prev_data + pf_data + [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)]

dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Moderately'
dm.params['global_priors']['smoothness']['incidence']['amount'] = 'No Prior'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)




dm = initialize_model()
dm.description = 'increasing and decreasing priors to keep edges looking ok'
prev_data = [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
for d in prev_data:
    d['effective_sample_size'] = 1.e5
dm.params['delta_effect_prevalence'] = dict(mean=5.,
                                            std=.025)
dm.params['global_priors']['decreasing']['prevalence'] = dict(age_start=70, age_end=100)
dm.params['global_priors']['increasing']['excess_mortality'] = dict(age_start=70, age_end=100)

pf_data = [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)]
for d in pf_data + prev_data:
    d['effective_sample_size'] = 1.e5

dm.data = prev_data + pf_data + [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)]

dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Moderately'
dm.params['global_priors']['smoothness']['incidence']['amount'] = 'No Prior'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=True, store_results=False)
models.append(dm)




dm = initialize_model()
dm.description = 'best version before getting updated data'
prev_data = [d for d in dm.data if d.get('country_iso3_code') == 'GBR']
for d in prev_data:
    d['effective_sample_size'] = 1.e5
dm.params['delta_effect_prevalence'] = dict(mean=3.,
                                            std=.25)
dm.params['global_priors']['decreasing']['prevalence'] = dict(age_start=70, age_end=100)

pf_data = [d for d in dm.data if dm.relevant_to(d, 'prevalence_x_excess-mortality', region, year, sex)]
for d in pf_data + prev_data:
    d['effective_sample_size'] = 1.e5

dm.data = prev_data + pf_data + [d for d in dm.data if dm.relevant_to(d, 'all-cause_mortality', region, year, sex)]

dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Moderately'
dm.params['global_priors']['smoothness']['incidence']['amount'] = 'Moderately'
fit_posterior.fit_posterior(dm, region, sex, year, map_only=False, store_results=False)
models.append(dm)



for ii, dm in enumerate(models):
    book_graphics.plot_age_patterns(dm, region, year, sex, rate_types='prevalence excess-mortality prevalence_x_excess-mortality'.split(),
                                    yticks={'prevalence': [0, .0006, .0012], 'excess-mortality': [0, 8, 16], 'prevalence_x_excess-mortality': [0, .001, .002]})
    pl.figtext(0, .005,'(%s) %s' % ('abcdefghijklmnopqrstuvqxyz'[ii], dm.description), fontsize=10)
    pl.savefig('cirrhosis-%d.png'%ii)
    pl.savefig('cirrhosis-%d.pdf'%ii)

pl.show()
