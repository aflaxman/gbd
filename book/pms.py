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

region='world'
year=2005
sex='female'
keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])

print 'initializing model vars... ',
dm.calc_effective_sample_size(dm.data)
for k in keys:
    dm.fit_initial_estimate(k)

### @export 'empirical-priors'
dm.vars = dismod3.generic_disease_model.setup(dm, '%s+world+2005+female')

verbose=1
def map_fit(stoch_names):
    print '\nfitting', ' '.join(stoch_names)
    key = dismod3.utils.gbd_key_for('%s', region, year, sex)

    # fit all parameters besides over-dispersion
    map = mc.MAP([dm.vars[key%type][subkey] for type in stoch_names for subkey in dm.vars[key%type] if not subkey.startswith('log_dispersion')])
    try:
        map.fit(method='fmin_powell', verbose=verbose)
    except KeyboardInterrupt:
        print 'User halted optimization routine before optimal value found'

    for type in stoch_names:
        for subkey in ['age_coeffs_mesh', 'dispersion']:
            if subkey in dm.vars[key%type]:
                print key%type, subkey, pl.atleast_1d(dm.vars[key%type][subkey].value).round(0)

    return map

# use map as initial values
print 'initializing MAP object... ',
dm.map = map_fit('incidence remission excess-mortality mortality relative-risk smr prevalence_x_excess-mortality duration bins prevalence'.split())
print 'initialization completed'

# fit with mcmc
dm.mcmc = mc.MCMC(dm.vars)
for k in keys:
    if 'dispersion_step_sd' in dm.vars[k]:
        dm.mcmc.use_step_method(mc.Metropolis, dm.vars[k]['log_dispersion'],
                                proposal_sd=dm.vars[k]['dispersion_step_sd'])
    if 'age_coeffs_mesh_step_cov' in dm.vars[k]:
        dm.mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['age_coeffs_mesh'],
                                cov=dm.vars[k]['age_coeffs_mesh_step_cov'], verbose=0)
dm.mcmc.sample(iter=5000, burn=3000, thin=10, verbose=verbose)

### @export 'save'
book_graphics.save_json('pms.json', vars())
