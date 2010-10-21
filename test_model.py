""" Script to test the integrated systems model for disease

"""

from pylab import *
from dismod3.disease_json import DiseaseJson
from dismod3 import neg_binom_model
import dismod3.utils

def test_single_rate():
    """ Test fit for a single low-noise data point"""

    # load model to test fitting
    dm = DiseaseJson(file('tests/single_low_noise.json').read())

    # fit empirical priors
    from dismod3 import neg_binom_model
    neg_binom_model.fit_emp_prior(dm, 'prevalence', '/dev/null')

    # compare fit to data
    check_emp_prior_fits(dm)

    return dm
    assert abs(data_prediction / dm.value_per_1(d) - 1.) < .1, 'Prediction should be closer to data'


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
    gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1)     ## then sample the posterior via MCMC

    # compare fit to data
    check_posterior_fits(dm)


def test_dismoditis():
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

    # fit posterior
    delattr(dm, 'vars')  # remove vars so that gbd_disease_model creates its own version
    from dismod3 import gbd_disease_model
    keys = dismod3.utils.gbd_keys(region_list=['asia_southeast'],
                                  year_list=[1990], sex_list=['male'])
    gbd_disease_model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1)     ## then sample the posterior via MCMC

    # compare fit to data
    check_posterior_fits(dm)


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
    gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1)     ## then sample the posterior via MCMC

    # compare fit to data
    check_posterior_fits(dm)

def check_posterior_fits(dm):
    print '*********************'
    for d in dm.data:
        type = d['parameter'].split()[0]
            
        prediction = dm.get_mcmc('mean', dismod3.utils.gbd_key_for(type, d['gbd_region'], d['year_start'], d['sex']))

        data_prediction = dismod3.utils.rate_for_range(prediction,
                                                       arange(d['age_start'], d['age_end']+1),
                                                       d['age_weights'])

        # test distance of predicted data value from observed data value
        print type, d['age_start'], dm.value_per_1(d), data_prediction, abs(100 * (data_prediction / dm.value_per_1(d) - 1.))
        #assert abs((.01 + data_prediction) / (.01 + dm.value_per_1(d)) - 1.) < 1., 'Prediction should be closer to data'
    print '*********************\n\n\n\n\n'


def check_emp_prior_fits(dm):
    # compare fit to data
    print '*********************'
    for d in dm.vars['data']:
        type = d['parameter'].split()[0]
        prior = dm.get_empirical_prior(type)     
        prediction = neg_binom_model.predict_country_rate(dismod3.utils.gbd_key_for(type, 'asia_southeast', 1990, 'male'), d['country_iso3_code'],
                                                          prior['alpha'], prior['beta'], prior['gamma'], dm.get_covariates(), lambda f, age: f, arange(101))
        data_prediction = dismod3.utils.rate_for_range(prediction,
                                                       arange(d['age_start'], d['age_end']+1),
                                                       d['age_weights'])

        # test distance of predicted data value from observed data value
        print type, d['age_start'], dm.value_per_1(d), data_prediction, abs(100 * (data_prediction / dm.value_per_1(d) - 1.))
        #assert abs((.01 + data_prediction) / (.01 + dm.value_per_1(d)) - 1.) < 1., 'Prediction should be closer to data'
    print '*********************\n\n\n\n\n'


if __name__ == '__main__':
    test_dismoditis_wo_prevalence()
    test_triangle_pattern()
    test_linear_pattern()
    test_single_rate()
    test_dismoditis()
