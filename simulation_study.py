""" Script to test the integrated systems model for disease

"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 


from pylab import *
import pymc as mc
import random    

from dismod3.disease_json import DiseaseJson
from dismod3 import neg_binom_model
import dismod3.utils
from test_model import summarize_acorr, check_emp_prior_fits, check_posterior_fits


def fit_simulated_disease(n=300, cv=2.):
    """ Test fit for simulated disease data with noise and missingness"""

    # load model to test fitting
    dm = DiseaseJson(file('tests/simulation_gold_standard.json').read())
    
    # adjust any priors and covariates as desired
    dm.set_param_age_mesh(arange(0,101,2))
    for type in 'incidence prevalence remission excess_mortality'.split():
        dm.params['global_priors']['heterogeneity'][type] = 'Very'
    
    # filter and noise up data
    mort_data = []
    all_data = []
    for d in dm.data:
        d['truth'] = d['value']
        d['age_weights'] = array([1.])
        if d['data_type'] == 'all-cause mortality data':
            mort_data.append(d)
        else:
            if d['value'] > 0:
                se = (cv / 100.) * d['value']
                d['value'] = mc.rtruncnorm(d['truth'], se**-2, 0, np.inf)
                d['standard_error'] = se

            all_data.append(d)
    dm.data = random.sample(all_data, n) + mort_data

    # fit empirical priors and compare fit to data
    from dismod3 import neg_binom_model
    for rate_type in 'prevalence incidence remission excess-mortality'.split():
        neg_binom_model.fit_emp_prior(dm, rate_type, iter=10000, thin=5, burn=5000, dbname='/dev/null')
        check_emp_prior_fits(dm)


    # fit posterior
    delattr(dm, 'vars')  # remove vars so that gbd_disease_model creates its own version
    from dismod3 import gbd_disease_model
    keys = dismod3.utils.gbd_keys(region_list=['north_america_high_income'],
                                  year_list=[1990], sex_list=['male'])
    gbd_disease_model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=10000, thin=5, burn=5000, verbose=1, dbname='/dev/null')     ## then sample the posterior via MCMC


    print 'error compared to the noisy data (coefficient of variation = %.2f)' % cv
    check_posterior_fits(dm)

    dm.data = all_data
    for d in dm.data:
        d['value'] = d['truth']

    print 'error compared to the truth'
    are, coverage = check_posterior_fits(dm)
    print
    print 'Median Absolute Relative Error of Posterior Predictions:', median(are)
    print 'Pct coverage:', 100*mean(coverage)
    f = open('score_%d_%f.txt' % (n, cv), 'a')
    f.write('%10.10f,%10.10f\n' % (median(are), mean(coverage)))
    f.close()
    
    # generate plots of results
    region = 'north_america_high_income'
    year = 1990
    sex = 'male'
    keys = dismod3.utils.gbd_keys(region_list=['north_america_high_income'], year_list=[1990], sex_list=['male'])
    dismod3.tile_plot_disease_model(dm, keys, defaults={})
    savefig('sim_%d_%f-posterior.png' % (n, cv))
    for k in dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex]):
        if dm.vars[k].get('data'):
            dismod3.plotting.plot_posterior_predicted_checks(dm, k)
            savefig('sim_%d_%f-check.png' % (n, cv))

    return dm


if __name__ == '__main__':
    import optparse
    usage = 'usage: %prog [options] n cv'
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()

    if len(args) != 2:
        # print summary results
        print '%25s\t|\t%s\t|\t%s\t|\t%s\t\\\\' % ('fname', 'n', 'MARE (%)', 'Coverage (%)')
        import glob
        for fname in sorted(glob.glob('score*.txt')):
            X = csv2rec(fname, names=['mare', 'coverage'])
            print '%25s\t|\t%d\t|\t%2.2f\t|\t%.2f\t\\\\' % (fname, len(X), mean(X.mare), mean(X.coverage)*100)
    else:

        try:
            n = int(args[0])
        except ValueError:
            parser.error('n must be an integer')

        try:
            cv = float(args[1])
            assert cv > 0, 'cv must be positive'
        except ValueError:
            parser.error('cv must be an integer')


        dm = fit_simulated_disease(n, cv)
