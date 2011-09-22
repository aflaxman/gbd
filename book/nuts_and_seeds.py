import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

type = 'remission'
region = 'europe_western'
year = '1990'
sex = 'male'
ymax = .002

def initialize_model():
    ### @export 'load model'
    dm = dismod3.load_disease_model(19271)
    ### @export 'initialize model data'
    dm.data = [d for d in dm.data if dm.relevant_to(d, type, region, year, sex)]
    for d in dm.data:
        d['standard_error'] = float(d['sd_1enadj'] or d['parameter_value_old'])/10000. / pl.sqrt(d['effective_sample_size'])
        d.pop('effective_sample_size')
    # fit model
    dm.clear_fit()
    dm.clear_empirical_prior()
    dismod3.neg_binom_model.covariate_hash = {}
    return dm

models = []

dm = initialize_model()
dm.params['gamma_effect_%s'%type] = dict(mean=list(pl.log(.00001*pl.ones_like(dm.get_estimate_age_mesh()))),
                                         std=list(2.*pl.ones_like(dm.get_estimate_age_mesh())))
dismod3.neg_binom_model.fit_emp_prior(dm, type, map_only=True, store_results=False)
models.append(dm)

dm = initialize_model()
dm.params['gamma_effect_%s'%type] = dict(mean=list(pl.log(.0001*pl.ones_like(dm.get_estimate_age_mesh()))),
                                         std=list(2.*pl.ones_like(dm.get_estimate_age_mesh())))
dismod3.neg_binom_model.fit_emp_prior(dm, type, map_only=True, store_results=False)
models.append(dm)

dm = initialize_model()
dm.params['gamma_effect_%s'%type] = dict(mean=list(pl.log(.001*pl.ones_like(dm.get_estimate_age_mesh()))),
                                         std=list(2.*pl.ones_like(dm.get_estimate_age_mesh())))
dismod3.neg_binom_model.fit_emp_prior(dm, type, map_only=True, store_results=False)
models.append(dm)

dm = initialize_model()
dm.params['gamma_effect_%s'%type] = dict(mean=list(pl.log(.01*pl.ones_like(dm.get_estimate_age_mesh()))),
                                         std=list(2.*pl.ones_like(dm.get_estimate_age_mesh())))
dismod3.neg_binom_model.fit_emp_prior(dm, type, map_only=True, store_results=False)
models.append(dm)


dm = initialize_model()
dismod3.neg_binom_model.fit_emp_prior(dm, type, map_only=True, store_results=False)
models.append(dm)

for ii, dm in enumerate(models):
    pl.subplot(1, len(models), len(models)-ii)
    dismod3.plotting.plot_intervals(dm, [d for d in dm.data if dm.relevant_to(d, type, region, year, sex)],
                                    color='black', print_sample_size=False, alpha=1., plot_error_bars=False,
                                    linewidth=2)
    key = dismod3.utils.gbd_key_for(type, region, year, sex)
    dm.vars = {key: dm.vars} # HACK: put the rate model in the dictionary as expected
    book_graphics.plot_rate(dm, key)
    dm.vars = dm.vars[key]
    pl.axis([0, 100, 0, ymax])
    pl.title('')
    pl.xticks([0,25,50,75])
    pl.xlabel('Age (Years)')
    pl.text(5, ymax*.9, '$\\gamma\sim N(%.2f, 2^2)$' % dm.params.get('gamma_effect_%s'%type, {}).get('mean', [0])[0], va='top', ha='left')
        
pl.title('Nuts and Seeds with different priors on $\\gamma$')
pl.show()
