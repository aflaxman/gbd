""" Script to fit the North America High Income Hep C data for 1990

"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

import pandas
import pylab as pl
import pymc as mc
    
from dismod3.disease_json import DiseaseJson
import dismod3.utils

import data
import data_model
import fit_model
import covariate_model
import graphics

reload(data_model)
reload(graphics)

pl.seterr('ignore')

id = 8788
#id = 23884

def summarize(name, df):
    bias = (df['value'] - df['mu_pred']).mean()
    mae = pl.median(pl.absolute((df['value'] - df['mu_pred'])))
    pc = pl.mean(pl.absolute(df['value'] - df['mu_pred']) < 1.96*df['sigma_pred'])
    return pandas.DataFrame(dict(name=[name], bias=[bias], mae=[mae], pc=[pc*100]))

def hold_out_quality():
    import glob

    results = pandas.DataFrame()
    for fname in sorted(glob.glob('hep_c_figs/*.csv')):
        df = pandas.read_csv(fname)
        results = results.append(summarize(fname, df[df['effective_sample_size']==0]), ignore_index=True)
    results = results.reindex(columns='name bias mae pc'.split())
    results = results.sort('mae')
    print results
    

def store_fit(dm, key, est_k):
    pl.figure()
    graphics.plot_data_bars(dm.model.input_data)

    if key.startswith('sim_'):
        dismod3.plotting.plot_fit(dm, 'mcmc_mean', key.replace('sim_', ''), color='blue', linewidth=3, label='Truth')
        pl.plot([0], [0], color='red', linewidth=3, label='Estimate')
    else:
        dismod3.plotting.plot_mcmc_fit(dm, key, color='blue')
        if id == 23884:
            pl.plot([0], [0], color='blue', linewidth=3, label='Standard Hierarchy')
            pl.plot([0], [0], color='red', linewidth=3, label='Custom Hierarchy')
        elif id == 8788:
            pl.plot([0], [0], color='blue', linewidth=3, label='Old')
            pl.plot([0], [0], color='red', linewidth=3, label='New')

    est_k.sort(axis=0)
    dm.set_mcmc('mean', key, pl.mean(est_k, axis=0))
    dm.set_mcmc('median', key, pl.median(est_k, axis=0))
    n = len(est_k)
    dm.set_mcmc('lower_ui', key, est_k[.025*n,:])
    dm.set_mcmc('upper_ui', key, est_k[.975*n,:])

    dismod3.plotting.plot_mcmc_fit(dm, key, color='red')
    pl.title(key.replace('+', ', '))

    #graphics.plot_one_ppc(dm.vars, key.replace('+', ', '))
    #graphics.plot_convergence_diag(dm.vars)
    pl.legend(shadow=True, fancybox=True)
    graphics.expand_axis(.05)
    pl.savefig('hep_c_figs/%s.png'%key)

def hep_c_fit(regions, prediction_years, data_year_start=-pl.inf, data_year_end=pl.inf, egypt_flag=False):
    """ Fit prevalence for regions and years specified """
    print '\n***************************\nfitting %s for %s (using data from years %f to %f)' % (regions, prediction_years, data_year_start, data_year_end)
    
    ## load model to fit
    dm = dismod3.load_disease_model(id)
    dm.data = [d for d in dm.data if d['data_type'] == 'prevalence data']

    ## adjust the expert priors
    dm.params['global_priors']['heterogeneity']['prevalence'] = 'Moderately'
    dm.params['global_priors']['smoothness']['prevalence']['amount'] = 'Moderately'
    dm.params['global_priors']['level_value']['prevalence']['age_before'] = 1
    dm.params['global_priors']['decreasing']['prevalence'] = dict(age_start=55, age_end=100)
    dm.params['global_priors']['parameter_age_mesh'] = [0, 1, 15, 25, 35, 45, 55, 100]

    # include a study-level covariate for 'bias'
    covariates_dict = dm.get_covariates()
    covariates_dict['Study_level']['bias']['rate']['value'] = 0
    covariates_dict['Study_level']['bias']['error']['value'] = 1
    covariates_dict['Country_level'] = {}

    ## select relevant prevalence data
    # TODO: streamline data selection functions
    if egypt_flag:
        dm.data = [d for d in dm.data if d['country_iso3_code'] == 'EGY']
    else:
        dm.data = [d for d in dm.data if
                   dismod3.utils.clean(d['gbd_region']) in regions
                   and float(d['year_end']) >= data_year_start
                   and float(d['year_start']) <= data_year_end
                   and d['country_iso3_code'] != 'EGY']

    ## create, fit, and save rate model
    import simplejson as json
    dm.model = data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))

    # uncomment following lines to hold out random 25% of observations for cross-validation
    #import random
    #i = random.sample(dm.model.input_data.index, len(dm.model.input_data)/4)
    #dm.model.input_data['effective_sample_size'][i] = 0.

    # add rows to the output template for sex, year == total, all
    total_template = dm.model.output_template.groupby('area').mean()
    total_template['area'] = total_template.index
    total_template['sex'] = 'total'
    total_template['year'] = 'all'
    dm.model.output_template = dm.model.output_template.append(total_template, ignore_index=True)

    dm.vars = data_model.data_model('p', dm.model, 'p',
                                    root_area='all', root_sex='total', root_year='all',
                                    mu_age=None,
                                    mu_age_parent=None,
                                    sigma_age_parent=None)

    #dm.map, dm.mcmc = fit_model.fit_data_model(dm.vars, 105, 0, 1, 100)
    dm.map, dm.mcmc = fit_model.fit_data_model(dm.vars, 4040, 2000, 20, 100)

    # add prediction values to DataFrame for posterior predictive check
    dm.vars['data']['mu_pred'] = dm.vars['p_pred'].stats()['mean']
    dm.vars['data']['sigma_pred'] = dm.vars['p_pred'].stats()['standard deviation']
    dm.vars['data']['residual'] = dm.vars['data']['value'] - dm.vars['data']['mu_pred']
    dm.vars['data']['abs_residual'] = pl.absolute(dm.vars['data']['residual'])
    dm.vars['data']['logp'] = [mc.negative_binomial_like(n*p_obs, n*p_pred+1.e-3, (n*p_pred+1.e-3)*d) for n, p_obs, p_pred, d \
                                   in zip(dm.vars['data']['effective_sample_size'], dm.vars['data']['value'], dm.vars['data']['mu_pred'], dm.vars['delta'].stats()['mean'])]

    dm.vars['data'] = dm.vars['data'].sort('logp')
    print dm.vars['data'].filter('area sex year_start age_start age_end effective_sample_size value mu_pred logp'.split())
    dm.vars['data'].filter('area sex year_start age_start age_end effective_sample_size value mu_pred sigma_pred logp'.split()).to_csv('hep_c_figs/%s.csv'%'+'.join([str(x) for x in regions + prediction_years]))


    keys = dismod3.utils.gbd_keys(type_list=['prevalence'],
                                  region_list=regions,
                                  year_list=prediction_years)

    for key in keys:
        t, r, y, s = dismod3.utils.type_region_year_sex_from_key(key)
        # special case, make sure EGY is capitalized
        if r == 'egy':
            r = 'EGY'
        est_k = covariate_model.predict_for(dm.model, 'all', 'total', 'all', r, s, int(y), 1., dm.vars, 0., 1.)
        store_fit(dm, key, est_k)
    print keys
    sim_and_fit(dm, keys)
 
    return dm

def sim_and_fit(dm, keys):
    data = dm.vars['data'].copy()

    delta_true = .1
    p = data['mu_pred']+1.e-6
    n = data['effective_sample_size']*10 # scale up ess, to simulate studies that have small sampling error (expect to see less zeros in data than unscaled)

    data['true'] = p
    data['value'] = (1.0 * mc.rnegative_binomial(n*p, delta_true*n*p) )/ n

    dm.model.input_data = data
    
    dm.sim_vars = data_model.data_model('p', dm.model, 'p',
                                        root_area='all', root_sex='total', root_year='all',
                                        mu_age=None,
                                        mu_age_parent=None,
                                        sigma_age_parent=None)

    print dm.sim_vars['data'].filter('area sex year_start age_start age_end effective_sample_size true value'.split())

    #dm.sim_map, dm.sim_mcmc = fit_model.fit_data_model(dm.sim_vars, 105, 0, 1, 100)
    dm.sim_map, dm.sim_mcmc = fit_model.fit_data_model(dm.sim_vars, 4040, 2000, 20, 100)

    for key in keys:
        t, r, y, s = dismod3.utils.type_region_year_sex_from_key(key)
        # special case, make sure EGY is capitalized
        if r == 'egy':
            r = 'EGY'
        est_k = covariate_model.predict_for(dm.model, 'all', 'total', 'all', r, s, int(y), 1., dm.sim_vars, 0., 1.)

        results = pandas.DataFrame(dict(true=dm.get_mcmc('mean', key), mu_pred=pl.mean(est_k, axis=0), sigma_pred=pl.std(est_k, axis=0)))
        results['residual'] = results['true'] - results['mu_pred']
        print '\nbias', results['residual'].mean(),
        print 'mae', pl.median(pl.absolute(results['residual'])),
        print 'pc', pl.mean(pl.absolute(results['residual']) < 1.96*results['sigma_pred'])

        results.to_csv('hep_c_figs/sim_%s.csv'%key)

        store_fit(dm, 'sim_' + key, est_k)


if __name__ == '__main__':

    dm = hep_c_fit('caribbean latin_america_tropical latin_america_andean latin_america_central latin_america_southern'.split(), [1990, 2005])
    dm = hep_c_fit('sub-saharan_africa_central sub-saharan_africa_southern sub-saharan_africa_west'.split(), [1990, 2005])
    
    for r in 'europe_eastern europe_central asia_central asia_east asia_south asia_southeast australasia oceania sub-saharan_africa_east asia_pacific_high_income'.split():
        dm = hep_c_fit([r], [1990, 2005])

    for r in 'europe_western north_america_high_income'.split():
        dm = hep_c_fit([r], [1990], data_year_end=1997)
        dm = hep_c_fit([r], [2005], data_year_start=1997)

    dm_egypt = hep_c_fit(['EGY'], [1990, 2005], egypt_flag=True)
    dm_na_me = hep_c_fit(['north_africa_middle_east'], [1990, 2005])

    # combine prevalence curves for egypt and rest of north africa
    dm = dismod3.load_disease_model(id)
    dm.data = [d for d in dm.data if d['data_type'] == 'prevalence data']
    dm.data = [d for d in dm.data if dismod3.utils.clean(d['gbd_region']) == 'north_africa_middle_east']
    covariates_dict = dm.get_covariates()
    covariates_dict['Country_level'] = {}
    import simplejson as json
    dm.model = data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))

    for y in [1990, 2005]:
        for s in ['male', 'female']:
            key = 'prevalence+EGY+%d+%s' % (y, s)
            from dismod3 import neg_binom_model
            prev_1 = covariate_model.predict_for(dm_egypt.model, 'EGY', 'total', 'all', 'EGY', s, y, 1., dm_egypt.vars, 0., 1.)
            pop_1 = neg_binom_model.population_by_age[('EGY', str(y), s)]

            key = 'prevalence+north_africa_middle_east+%d+%s' % (y, s)
            prev_0 = covariate_model.predict_for(dm_na_me.model, 'all', 'total', 'all', 'north_africa_middle_east', s, y, 1., dm_na_me.vars, 0., 1.)
            pop_0 = neg_binom_model.regional_population(key)

            # generate population weighted average
            est_k = (prev_0 * (pop_0 - pop_1) + prev_1 * pop_1) / pop_0
            store_fit(dm, key, est_k)

            # save results
            #dismod3.post_disease_model(dm_na_me)
