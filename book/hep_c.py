import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}

regions = ['north_america_high_income']
#regions = ['europe_western']
data_year_start = 1980; data_year_end = 1997; prediction_years = [1990]
data_year_start = 1997; data_year_end = 2011; prediction_years = [2005]

### @export 'data'
# load model
dm = dismod3.load_disease_model(16391)

# set expert priors and other model parameters
dm.params['global_priors']['heterogeneity']['prevalence'] = 'Slightly'
dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'No Prior'
dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Slightly'
#dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Moderately'
#dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Very'
dm.params['global_priors']['level_value']['prevalence']['age_before'] = 1

dm.params['covariates']['Study_level']['bias']['rate']['value'] = 1
dm.params['covariates']['Study_level']['bias']['error']['value'] = 1  # TODO: ensure that this is doing what I think it is doing.
for cv in dm.params['covariates']['Country_level']:
    dm.params['covariates']['Country_level'][cv]['rate']['value'] = 0

# clear any fit and priors
dm.clear_fit()
dm.clear_empirical_prior()
dismod3.neg_binom_model.covariate_hash = {}

# initialize model data
dm.data = [d for d in dm.data if
           d['data_type'] == 'prevalence data'
           and dismod3.utils.clean(d['gbd_region']) in regions
           and float(d['year_end']) >= data_year_start
           and float(d['year_start']) <= data_year_end
           and d['country_iso3_code'] != 'EGY']

dm.vars = {}
keys = dismod3.utils.gbd_keys(type_list=['prevalence'],
                              region_list=regions,
                              year_list=prediction_years)
k = keys[0]  # looks like k='prevalence+asia_south+1990+male'
dm.vars[k] = dismod3.neg_binom_model.setup(dm, k, dm.data)


# fit the model
def map_fit(stoch_names):
    print '\nfitting', ' '.join(stoch_names)
    map = mc.MAP([dm.vars[k][key] for key in stoch_names] + [dm.vars[k]['observed_counts'], dm.vars[k]['rate_potential'], dm.vars[k]['priors']])
    try:
        map.fit(method='fmin_powell', verbose=verbose)
    except KeyboardInterrupt:
        debug('User halted optimization routine before optimal value found')
    for key in stoch_names:
        print key, dm.vars[k][key].value.round(2)
    sys.stdout.flush()

verbose = 1
stoch_names = 'region_coeffs age_coeffs_mesh study_coeffs'.split()
## start by optimizing parameters separately
for key in stoch_names:
    map_fit([key])
## then fit them all together
map_fit(stoch_names)
# now find the over-dispersion parameter that matches these values
map_fit(['log_dispersion'])


mcmc = mc.MCMC(dm.vars)


if 'dispersion_step_sd' in dm.vars[k]:
    mcmc.use_step_method(mc.Metropolis, dm.vars[k]['log_dispersion'],
                            proposal_sd=dm.vars[k]['dispersion_step_sd'])
if 'age_coeffs_mesh_step_cov' in dm.vars[k]:
    mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['age_coeffs_mesh'],
                            cov=dm.vars[k]['age_coeffs_mesh_step_cov'], verbose=0)

    # TODO: make a wrapper function for handling this adaptive metropolis setup
    stoch_list = [dm.vars[k]['study_coeffs'], dm.vars[k]['region_coeffs'], dm.vars[k]['age_coeffs_mesh']]
    d1 = len(dm.vars[k]['study_coeffs'].value)
    d2 = len(dm.vars[k]['region_coeffs_step_cov'])
    d3 = len(dm.vars[k]['age_coeffs_mesh_step_cov'])
    C = pl.eye(d1+d2+d3)
    C[d1:(d1+d2), d1:(d1+d2)] = dm.vars[k]['region_coeffs_step_cov']
    C[(d1+d2):(d1+d2+d3), (d1+d2):(d1+d2+d3)] = dm.vars[k]['age_coeffs_mesh_step_cov']
    C *= .01
    mcmc.use_step_method(mc.AdaptiveMetropolis, stoch_list, cov=C)

    # more step methods
    mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['study_coeffs'])
    mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['region_coeffs'], cov=dm.vars[k]['region_coeffs_step_cov'])
    mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['age_coeffs_mesh'], cov=dm.vars[k]['age_coeffs_mesh_step_cov'])
mcmc.use_step_method(mc.AdaptiveMetropolis, stoch_list, cov=C)
mcmc.sample(iter=10000, burn=5000, thin=5, verbose=1)

k0 = k
for k in keys:
    # save the results in the disease model
    dm.vars[k] = dm.vars[k0]

    dismod3.neg_binom_model.store_mcmc_fit(dm, k, dm.vars[k])

    # generate plots of results
    dismod3.plotting.tile_plot_disease_model(dm, [k], defaults={'ymax':.15, 'alpha': .5})
    pl.savefig('hep_c-posterior-%s.png' % k)

# summarize fit quality graphically, as well as parameter posteriors
dismod3.plotting.plot_posterior_predicted_checks(dm, k0)

### @export 'save'
book_graphics.save_json('hep_c.json', vars())
