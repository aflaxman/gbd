""" Script to test the integrated systems model for disease

"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 


from pylab import *
import pymc as mc
import inspect
    
from dismod3.disease_json import DiseaseJson
from dismod3 import neg_binom_model
import dismod3.utils

def test_single_rate():
    """ Test fit for a single low-noise data point"""

    # load model to test fitting
    dm = DiseaseJson(file('tests/single_low_noise.json').read())

    # fit empirical priors
    neg_binom_model.fit_emp_prior(dm, 'prevalence', '/dev/null')

    # compare fit to data
    check_emp_prior_fits(dm)

def hep_c_fit():
    """ Test fit for prevalence only"""

    # load model to test fitting
    dm = DiseaseJson(file('tests/hep_c.json').read())

    # select only some prevalence data
    dm.data = [d for d in dm.data if dismod3.utils.clean(d['gbd_region']) == 'north_america_high_income']

    # fit empirical priors
    neg_binom_model.fit_emp_prior(dm, 'prevalence', '/dev/null')

    # fit posterior
    prior_dict = dm.get_empirical_prior('prevalence')


    # TODO: connect several prevalence vars with hierarchical potential
    keys = dismod3.utils.gbd_keys(type_list=['prevalence'],
                                  region_list=['north_america_high_income'],
                                  year_list=[1990],
                                  sex_list=['male'])
    from dismod3.gbd_disease_model import relevant_to
    from dismod3.utils import type_region_year_sex_from_key
    dm.vars = {}
    for k in keys:
        t,r,y,s = type_region_year_sex_from_key(k)

        data = [d for d in dm.data if relevant_to(d, 'all', r, y, s)]

        dm.vars[k] = neg_binom_model.setup(dm, k,
                                      data,
                                      emp_prior=prior_dict)
    
    # TODO: load initial values from empirical prior expectation

    import pymc as mc
    dbname = '/dev/null'
    dm.mcmc = mc.MCMC(dm.vars, db='pickle', dbname=dbname)
    for k in keys:
        if 'dispersion_step_sd' in dm.vars[k]:
            dm.mcmc.use_step_method(mc.Metropolis, dm.vars[k]['log_dispersion'],
                                    proposal_sd=dm.vars[k]['dispersion_step_sd'])
        if 'age_coeffs_mesh_step_cov' in dm.vars[k]:
            dm.mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['age_coeffs_mesh'],
                                    cov=dm.vars[k]['age_coeffs_mesh_step_cov'], verbose=0)

    iter = 1000
    thin = 5
    burn = iter*thin
    try:
        dm.mcmc.sample(iter=iter*thin+burn, thin=thin, burn=burn, verbose=1)
    except KeyboardInterrupt:
        # if user cancels with cntl-c, save current values for "warm-start"
        pass
    dm.mcmc.db.commit()


    # make map object to keep store function happy
    class my_obj:
        pass
    dm.map = my_obj()
    dm.map.AIC = np.nan
    dm.map.BIC = np.nan

    for k in keys:
        neg_binom_model.store_mcmc_fit(dm, k, dm.vars[k])

        # check autocorrelation to confirm chain has mixed
        summarize_acorr(dm.vars[k]['rate_stoch'].trace())

        # generate plots of results
        dismod3.tile_plot_disease_model(dm, [k], defaults={'ymax':.05})
        dm.savefig('dm-%d-posterior-%s.%f.png' % (dm.id, k, random()))

    return dm
                
def summarize_acorr(x):
    x = x - np.mean(x, axis=0)
    print '*********************', inspect.stack()[1][3]
    for a in np.arange(0,101,10):
        acorr5 = dot(x[5:, a], x[:-5, a]) / dot(x[5:, a], x[5:, a])
        acorr10 = dot(x[10:, a], x[:-10, a]) / dot(x[10:, a], x[10:, a])
        print 'a: %d, c5: %.2f, c10: %.2f' % (a, acorr5*100, acorr10*100)
    print '*********************'    

def test_simulated_disease():
    """ Test fit for simulated disease data"""

    # load model to test fitting
    dm = DiseaseJson(file('tests/test_disease_1.json').read())

    # filter and noise up data
    cov = .5
    
    data = []
    for d in dm.data:
        d['truth'] = d['value']
        if dismod3.utils.clean(d['gbd_region']) == 'north_america_high_income':
            if d['data_type'] == 'all-cause mortality data':
                data.append(d)
            else:
                se = (cov * d['value'])
                d['value'] = mc.rtruncnorm(d['truth'], se**-2, 0, np.inf)
                d['age_start'] -= 5
                d['age_end'] = d['age_start']+9
                d['age_weights'] = np.ones(d['age_end']-d['age_start']+1)
                d['age_weights'] /= float(len(d['age_weights']))

                d['standard_error'] = se

                data.append(d)

    dm.data = data
    
    # fit empirical priors and compare fit to data
    from dismod3 import neg_binom_model
    for rate_type in 'prevalence incidence remission excess-mortality'.split():
        neg_binom_model.fit_emp_prior(dm, rate_type, '/dev/null')
        check_emp_prior_fits(dm)


    # fit posterior
    delattr(dm, 'vars')  # remove vars so that gbd_disease_model creates its own version
    from dismod3 import gbd_disease_model
    keys = dismod3.utils.gbd_keys(region_list=['north_america_high_income'],
                                  year_list=[1990], sex_list=['male'])
    gbd_disease_model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1, dbname='/dev/null')     ## then sample the posterior via MCMC


    print 'error compared to the noisy data (coefficient of variation = %.2f)' % cov
    check_posterior_fits(dm)


    for d in dm.data:
        d['value'] = d['truth']
        d['age_start'] += 5
        d['age_end'] = d['age_start']
        d['age_weights'] = np.ones(d['age_end']-d['age_start']+1)
        d['age_weights'] /= float(len(d['age_weights']))

    print 'error compared to the truth'
    check_posterior_fits(dm)

    return dm

def test_mesh_refinement():
    """ Compare fit for coarse and fine age mesh"""

    # load model and fit it
    dm1 = DiseaseJson(file('tests/single_low_noise.json').read())
    dm1.set_param_age_mesh(arange(0,101,20))
    from dismod3 import neg_binom_model
    neg_binom_model.fit_emp_prior(dm1, 'prevalence', '/dev/null')

    # load another copy and fit it with a finer age mesh
    dm2 = DiseaseJson(file('tests/single_low_noise.json').read())
    dm2.set_param_age_mesh(arange(0,101,5))
    from dismod3 import neg_binom_model
    neg_binom_model.fit_emp_prior(dm2, 'prevalence', '/dev/null')

    # compare fits
    p1 = dm1.get_mcmc('emp_prior_mean', dismod3.utils.gbd_key_for('prevalence', 'asia_southeast', 1990, 'male'))
    p2 = dm2.get_mcmc('emp_prior_mean', dismod3.utils.gbd_key_for('prevalence', 'asia_southeast', 1990, 'male'))
    print p1[::20]
    print p2[::20]
    assert np.all(abs(p1[::20] / p2[::20] - 1.) < .05), 'Prediction should be closer to data'


def test_linear_pattern():
    """ Test fit for empirical prior to data showing a linearly increasing age pattern"""

    # load model to test fitting
    dm = DiseaseJson(file('tests/single_low_noise.json').read())

    # create linear age pattern data
    import copy
    d = dm.data.pop()
    for a in range(10, 100, 20):
        d = copy.copy(d)
        d['age_start'] = a
        d['age_end'] = a
        d['parameter_value'] = .01*a
        d['value'] = .01*a
        dm.data.append(d)

    # fit empirical priors
    from dismod3 import neg_binom_model
    neg_binom_model.fit_emp_prior(dm, 'prevalence', '/dev/null')

    # compare fit to data
    check_emp_prior_fits(dm)


def test_increasing_prior():
    """ Test fit for empirical prior to data showing a linearly increasing age pattern with a fine age mesh"""

    # load model to test fitting
    dm = DiseaseJson(file('tests/single_low_noise.json').read())

    dm.params['global_priors']['increasing']['incidence']['age_end'] = 100

    # create linear age pattern data
    import copy
    d = dm.data.pop()
    for a in range(10, 100, 10):
        d = copy.copy(d)
        d['age_start'] = a
        d['age_end'] = a
        d['parameter_value'] = .01*a
        d['value'] = .01*a
        dm.data.append(d)

    # fit empirical priors
    from dismod3 import neg_binom_model
    neg_binom_model.fit_emp_prior(dm, 'prevalence', '/dev/null')

    # compare fit to data, and check that it is increasing
    check_emp_prior_fits(dm)
    assert np.all(np.diff(dm.get_mcmc('emp_prior_mean', dismod3.utils.gbd_key_for('prevalence', 'asia_southeast', 1990, 'male'))) >= 0), 'expert prior says increasing'

def test_triangle_pattern():
    """ Test fit for empirical prior to data showing a linearly increasing age pattern"""

    # load model to test fitting
    dm = DiseaseJson(file('tests/single_low_noise.json').read())

    # create linear age pattern data
    import copy
    d = dm.data.pop()
    for a in range(10, 100, 20):
        d = copy.copy(d)
        d['age_start'] = a
        d['age_end'] = a
        d['parameter_value'] = .01*min(a, 100-a)
        d['value'] = .01*min(a, 100-a)
        dm.data.append(d)

    # fit empirical priors
    from dismod3 import neg_binom_model
    neg_binom_model.fit_emp_prior(dm, 'prevalence', '/dev/null')

    # compare fit to data
    check_emp_prior_fits(dm)

    # fit posterior
    delattr(dm, 'vars')  # remove vars so that gbd_disease_model creates its own version
    from dismod3 import gbd_disease_model
    keys = dismod3.utils.gbd_keys(region_list=['asia_southeast'],
                                  year_list=[1990], sex_list=['male'])
    gbd_disease_model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1, dbname='/dev/null')     ## then sample the posterior via MCMC

    # compare fit to data
    check_posterior_fits(dm)


def test_dismoditis():
    """ Test fit for simple example"""

    # load model to test fitting
    dm = DiseaseJson(file('tests/dismoditis.json').read())
    for d in dm.data:
        d['standard_error'] = .01
    # fit empirical priors
    neg_binom_model.fit_emp_prior(dm, 'prevalence', '/dev/null')
    check_emp_prior_fits(dm)
    neg_binom_model.fit_emp_prior(dm, 'incidence', '/dev/null')
    check_emp_prior_fits(dm)
    neg_binom_model.fit_emp_prior(dm, 'excess-mortality', '/dev/null')
    check_emp_prior_fits(dm)

    # fit posterior where there is no data
    delattr(dm, 'vars')  # remove vars so that gbd_disease_model creates its own version
    from dismod3 import gbd_disease_model
    keys = dismod3.utils.gbd_keys(region_list=['north_america_high_income'],
                                  year_list=[1990], sex_list=['male'])
    gbd_disease_model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1, dbname='/dev/null')     ## then sample the posterior via MCMC
    check_posterior_fits(dm)
    
    # check that prevalence is smooth near age zero
    prediction = dm.get_mcmc('mean', 'prevalence+north_america_high_income+1990+male')
    assert prediction[1]-prediction[0] < .01, 'prediction should be smooth near zero'

def test_hep_c():
    """ Test fit for subset of hep_c data

    data is filtered to include only prevalence with
    region == 'europe_western' and sex == 'all'
    """

    # load model to test fitting
    dm = DiseaseJson(file('tests/hep_c_europe_western.json').read())

    # fit empirical priors
    neg_binom_model.fit_emp_prior(dm, 'prevalence', '/dev/null')

    # fit posterior
    delattr(dm, 'vars')  # remove vars so that gbd_disease_model creates its own version
    from dismod3 import gbd_disease_model
    keys = dismod3.utils.gbd_keys(region_list=['europe_western'],
                                  year_list=[1990], sex_list=['male'])
    gbd_disease_model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1, dbname='/dev/null')     ## then sample the posterior via MCMC

    # check that prevalence is smooth near age zero
    prediction = dm.get_mcmc('mean', 'prevalence+europe_western+1990+male')
    print prediction
    assert prediction[100] < .1, 'prediction should not shoot up in oldest ages'

def test_opi():
    """ Test fit for subset of opi_dep data

    data is filtered to include only data for
    region == 'europe_central' and sex == 'male'
    """

    # load model to test fitting
    dm = DiseaseJson(file('tests/opi.json').read())

    dm.params['global_priors']['decreasing']['prevalence']['age_start'] = 60
    dm.params['global_priors']['decreasing']['prevalence']['age_end'] = 100

    fit_model(dm, 'europe_central', 1990, 'male')

    # check that prevalence is smooth near age zero
    prediction = dm.get_mcmc('mean', 'prevalence+europe_central+1990+male')
    print prediction
    assert prediction[80] > prediction[100], 'prediction should decrease at oldest ages'

    return dm

def test_ihd():
    """ Test fit for subset of ihd data

    data is filtered to include only data for
    region == 'europe_western' and sex == 'male'
    """

    # load model to test fitting
    dm = DiseaseJson(file('tests/ihd.json').read())

    fit_model(dm, 'europe_western', 1990, 'male')

    # check that prevalence is smooth around age 90
    prediction = dm.get_mcmc('mean', 'prevalence+europe_western+1990+male')
    print prediction
    assert prediction[89]/prediction[90] < .05, 'prediction should not change greatly at age 90'

    return dm

def fit_model(dm, region, year, sex):
    """ Fit the empirical priors, and the posterior for a specific region/year/sex
    """
    
    # fit empirical priors
    for rate_type in 'prevalence incidence remission excess-mortality'.split():
        neg_binom_model.fit_emp_prior(dm, rate_type, '/dev/null')

    # fit posterior
    delattr(dm, 'vars')  # remove vars so that gbd_disease_model creates its own version
    from dismod3 import gbd_disease_model
    keys = dismod3.utils.gbd_keys(region_list=[region],
                                  year_list=[year], sex_list=[sex])
    gbd_disease_model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1, dbname='/dev/null')     ## then sample the posterior via MCMC



def test_dismoditis_w_high_quality_data():
    """ Test fit for simple example"""

    # load model to test fitting
    dm = DiseaseJson(file('tests/dismoditis.json').read())

    # fit empirical priors
    neg_binom_model.fit_emp_prior(dm, 'prevalence', '/dev/null')
    check_emp_prior_fits(dm)
    neg_binom_model.fit_emp_prior(dm, 'incidence', '/dev/null')
    check_emp_prior_fits(dm)
    neg_binom_model.fit_emp_prior(dm, 'excess-mortality', '/dev/null')
    check_emp_prior_fits(dm)

    # fit posterior where there is data
    delattr(dm, 'vars')  # remove vars so that gbd_disease_model creates its own version
    from dismod3 import gbd_disease_model
    keys = dismod3.utils.gbd_keys(region_list=['asia_southeast'],
                                  year_list=[1990], sex_list=['male'])
    gbd_disease_model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1, dbname='/dev/null')     ## then sample the posterior via MCMC

    # compare fit to data
    check_posterior_fits(dm)

    return dm
    


def test_dismoditis_wo_prevalence():
    """ Test fit for simple example"""

    # load model to test fitting
    dm = DiseaseJson(file('tests/dismoditis.json').read())

    # remove all prevalence data
    dm.data = [d for d in dm.data if d['parameter'] != 'prevalence data']

    # fit empirical priors
    neg_binom_model.fit_emp_prior(dm, 'incidence', '/dev/null')
    check_emp_prior_fits(dm)
    neg_binom_model.fit_emp_prior(dm, 'excess-mortality', '/dev/null')
    check_emp_prior_fits(dm)

    # fit posterior
    delattr(dm, 'vars')  # remove vars so that gbd_disease_model creates its own version
    from dismod3 import gbd_disease_model
    keys = dismod3.utils.gbd_keys(region_list=['asia_southeast'],
                                  year_list=[1990], sex_list=['male'])
    #gbd_disease_model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1, dbname='/dev/null')     ## then sample the posterior via MCMC

    # compare fit to data
    check_posterior_fits(dm)

def check_posterior_fits(dm):
    print '*********************', inspect.stack()[1][3]
    for d in dm.data:
        type = d['data_type'].replace(' data', '')
            
        prediction = dm.get_mcmc('mean', dismod3.utils.gbd_key_for(type, d['gbd_region'], d['year_start'], d['sex']))
        if len(prediction) == 0:
            continue

        data_prediction = dismod3.utils.rate_for_range(prediction,
                                                       arange(d['age_start'], d['age_end']+1),
                                                       d['age_weights'])
        
        # test distance of predicted data value from observed data value
        print type, d['age_start'], dm.value_per_1(d), data_prediction, abs(100 * (data_prediction / dm.value_per_1(d) - 1.))
        #assert abs((.01 + data_prediction) / (.01 + dm.value_per_1(d)) - 1.) < 1., 'Prediction should be closer to data'
    print '*********************\n\n\n\n\n'


def check_emp_prior_fits(dm):
    # compare fit to data
    print '*********************', inspect.stack()[1][3]
    for d in dm.vars['data']:
        type = d['data_type'].replace(' data', '')
        prior = dm.get_empirical_prior(type)     
        prediction = neg_binom_model.predict_country_rate(dismod3.utils.gbd_key_for(type, 'asia_southeast', 1990, 'male'), d['country_iso3_code'],
                                                          prior['alpha'], prior['beta'], prior['gamma'], dm.get_covariates(), lambda f, age: f, arange(101))
        data_prediction = dismod3.utils.rate_for_range(prediction,
                                                       arange(d['age_start'], d['age_end']+1),
                                                       d['age_weights'])

        # test distance of predicted data value from observed data value
        print type, d['age_start'], dm.value_per_1(d), data_prediction, abs(100 * (data_prediction / dm.value_per_1(d) - 1.))
        #assert abs((.001 + data_prediction) / (.001 + dm.value_per_1(d)) - 1.) < .05, 'Prediction should be closer to data'
    print '*********************\n\n\n\n\n'


if __name__ == '__main__':
    for test in [
        test_ihd,
        test_opi,
        test_hep_c,
        test_dismoditis,
        test_dismoditis_w_high_quality_data,
        test_mesh_refinement,
        test_increasing_prior,
        test_dismoditis_wo_prevalence,
        test_triangle_pattern,
        test_linear_pattern,
        test_single_rate,
        ]:
        try:
            test()
        except AssertionError, e:
            print 'TEST FAILED', test
            print e
