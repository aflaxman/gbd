import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}

### @export 'data'
dm = dismod3.load_disease_model(16314)
dm.clear_fit()
dm.clear_empirical_prior()

dm.set_param_age_mesh([0, 15, 20, 25, 30, 35, 40, 45, 47, 49, 51, 53, 55, 100])

dm.calc_effective_sample_size(dm.data)
some_data = ([d for d in dm.data
              if d['data_type'] == 'prevalence data'
              and d['sex'] == 'female'
              and d['effective_sample_size'] > 1])

# TODO: replace fake year data with real year data (called Year Start (original) and End)

countries = pl.unique([s['region'] for s in some_data])
min_year = min([s['year_start'] for s in some_data])
max_year = max([s['year_end'] for s in some_data])
cy = ['%s-%d'%(s['region'], s['year_start']) for s in some_data]

n = pl.array([s['effective_sample_size'] for s in some_data])
r = pl.array([dm.value_per_1(s) for s in some_data])
s = pl.sqrt(r * (1-r) / n)
min_rate_per_100 = '%d' % round(min([dm.value_per_1(d) for d in some_data])*100)

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

print 'initializing MAP object... ',
map_fit(['incidence', 'bins', 'prevalence'])
map_fit(['remission'])
map_fit(['excess-mortality'])
dm.map = map_fit('incidence remission excess-mortality mortality relative-risk smr prevalence_x_excess-mortality duration bins prevalence'.split())

print 'initialization completed'

dm.mcmc = mc.MCMC(dm.vars)
age_stochs = []
for k in keys:
    if 'dispersion_step_sd' in dm.vars[k]:
        dm.mcmc.use_step_method(mc.Metropolis, dm.vars[k]['log_dispersion'],
                                proposal_sd=dm.vars[k]['dispersion_step_sd'])
    if 'age_coeffs_mesh_step_cov' in dm.vars[k]:
        age_stochs.append(dm.vars[k]['age_coeffs_mesh'])
        dm.mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['age_coeffs_mesh'],
                                cov=dm.vars[k]['age_coeffs_mesh_step_cov'], verbose=0)
key = dismod3.utils.gbd_key_for('%s', region, year, sex)
try:
    dm.mcmc.sample(iter=20000, burn=10000, thin=10, verbose=verbose)
except KeyboardInterrupt:
    # if user cancels with cntl-c, save current values for "warm-start"
    pass

### @export 'save'
book_graphics.save_json('pms.json', vars())
