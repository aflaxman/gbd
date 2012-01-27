import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}

### @export 'data'
region = 'north_america_high_income'
sex = 'female'
year = '2005'

heterogeneity = ['Slightly', 'Very']

for ii in range(2):
    # load model
    dm = dismod3.load_disease_model(16370)

    # set expert priors and other model parameters
    dm.set_param_age_mesh([0,15,20,25,30,35,40,45,50,100])

    dm.params['global_priors']['level_value']['incidence']['age_before'] = 15
    dm.params['global_priors']['level_value']['incidence']['age_after'] = 50
    dm.params['global_priors']['smoothness']['incidence']['age_start'] = 15

    dm.params['global_priors']['level_value']['remission']['age_before'] = 40
    dm.params['global_priors']['level_bounds']['remission']['upper'] = 10.
    
    dm.params['global_priors']['level_value']['excess_mortality']['age_before'] = 101
    
    dm.params['global_priors']['level_value']['prevalence']['age_before'] = 15
    dm.params['global_priors']['level_value']['prevalence']['age_after'] = 50
    dm.params['global_priors']['smoothness']['prevalence']['age_start'] = 15
    dm.params['global_priors']['heterogeneity']['prevalence'] = heterogeneity[ii]
    
    dm.params['covariates']['Country_level']['LDI_id_Updated_7July2011']['rate']['value'] = 0
    # for cv in dm.params['covariates']['Study_level'].keys():
    #     if not cv.startswith('cv_'):
    #         use_flag = False
    #     else:
    #         use_flag = ii  # ii = 0 -> no covariates
        
    #     dm.params['covariates']['Study_level'][cv]['rate']['value'] = use_flag


    # clear any fit and priors
    dm.clear_fit()
    dm.clear_empirical_prior()
    dismod3.neg_binom_model.covariate_hash = {}

    import fit_posterior
    data = dm.data
    fit_posterior.fit_posterior(dm, region, sex, year, map_only=False, store_results=False)
    dm.data = data # put data back in
    results[ii] = dm


pl.figure(**book_graphics.quarter_page_params)
r = pl.array([dm.value_per_1(d) for d in data if dm.relevant_to(d, 'prevalence', region, year, sex)])
n = [d['effective_sample_size'] for d in data if dm.relevant_to(d, 'prevalence', region, year, sex)]
sorted_indices = r.argsort().argsort()

pl.errorbar(sorted_indices, r, yerr=1.96*pl.sqrt(r*(1-r)/n), fmt='ks', mew=1, ms=5, mec='white', label='Observed value')

n_samples=16
jitter = mc.rnormal(0, .025**-2, n_samples*len(r)).reshape(n_samples, len(r))
marker = ['k^', 'ko']
label = ['Posterior prediction without covariates', 'Posterior prediction with covariates']
# for i in sorted_indices:
#     for j in range(2):
#         vars = results[j].vars[dismod3.utils.gbd_key_for('prevalence', region, year, sex)]
#         try:
#             pl.plot(i+.25*(j+1)+jitter, vars['predicted_rates'].trace()[:n_samples, i], marker[j])

#         except Exception, e:
#             print e
#             pl.plot(i+.25*(j+1), vars['predicted_rates'].value[i], marker[j], label=label[j])
#             results['posterior_dispersion_%d'%j] = 'TK'

#for j in range(2):
#    pl.plot([-100], [0], marker[j], label=label[j]) # HACK: off-screen point with label to appear in legend

for j in range(2):
    vars = results[j].vars[dismod3.utils.gbd_key_for('prevalence', region, year, sex)]
    stats = vars['predicted_rates'].stats()
    y = stats['quantiles'][50]
    yerr = [y-stats['quantiles'][2.5], stats['quantiles'][97.5]-y]
    pl.errorbar(sorted_indices+.25*(j+1), y, yerr,
                mew=1, ms=5, mec='white',
                fmt=marker[j], label=label[j])
    
    results['posterior_dispersion_%d'%j] = vars['dispersion'].stats()
    results['posterior_beta_%d'%j] = vars['study_coeffs'].stats()

pl.xlabel('PMS prevalence data from systematic review')
pl.xticks([])
pl.yticks([0, .2, .4, .6, .8, 1])
pl.ylabel('Prevalence (per 1)')
pl.axis([-.5, 11.5,-.01,1.])
pl.legend(loc='upper left', fancybox=True, shadow=True, numpoints=1)

pl.savefig('pms_ppc.pdf')


book_graphics.save_json('pms_ppc.json', vars())
