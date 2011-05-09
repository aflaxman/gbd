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
        dm.params['covariates']['Country_level']['LDI_id']['rate']['value'] = 0
    
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
                Y_i = mc.rtruncnorm(d['truth'], se**-2, 0, np.inf)
                d['value'] = Y_i
                d['standard_error'] = se
                d['effective_sample_size'] = Y_i * (1-Y_i) / se**2


            all_data.append(d)
    sampled_data = random.sample(all_data, n) + mort_data
    dm.data = sampled_data

    # fit empirical priors and compare fit to data
    from dismod3 import neg_binom_model
    for rate_type in 'prevalence incidence remission excess-mortality'.split():
        #neg_binom_model.fit_emp_prior(dm, rate_type, iter=1000, thin=1, burn=0, dbname='/dev/null')
        neg_binom_model.fit_emp_prior(dm, rate_type, iter=30000, thin=15, burn=15000, dbname='/dev/null')
        check_emp_prior_fits(dm)


    # fit posterior
    delattr(dm, 'vars')  # remove vars so that gbd_disease_model creates its own version
    from dismod3 import gbd_disease_model
    keys = dismod3.utils.gbd_keys(region_list=['north_america_high_income'],
                                  year_list=[1990], sex_list=['male'])
    gbd_disease_model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=30000, thin=15, burn=15000, verbose=1, dbname='/dev/null')     ## then sample the posterior via MCMC
    #gbd_disease_model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=1, burn=0, verbose=1, dbname='/dev/null')     ## fast for dev


    print 'error compared to the noisy data (coefficient of variation = %.2f)' % cv
    check_posterior_fits(dm)

    dm.data = all_data
    for d in dm.data:
        if d['data_type'] != 'all-cause mortality data':
            d['noisy_data'] = d['value']
            d['value'] = d['truth']

    print 'error compared to the truth'
    are, coverage = check_posterior_fits(dm)
    print
    print 'Median Absolute Relative Error of Posterior Predictions:', median(are)
    print 'Pct coverage:', 100*mean(coverage)
    f = open('score_%d_%f.txt' % (n, cv), 'a')
    f.write('%10.10f,%10.10f\n' % (median(are), mean(coverage)))
    f.close()

    dm.all_data = all_data
    dm.data = sampled_data
    for d in dm.data:
        if d['data_type'] != 'all-cause mortality data':
            d['value'] = d['noisy_data']

    generate_figure(dm, n, cv)

    return dm

def generate_figure(dm, n, cv):
    figure(dpi=600, figsize=(8.5, 4))
    
    # generate plots of results
    region = 'north_america_high_income'
    year = 1990
    sex = 'male'

    reload(dismod3.plotting)
    all_data_hash = dismod3.plotting.GBDDataHash(dm.all_data)
    data_hash = dismod3.plotting.GBDDataHash(dm.data)

    label = { 'prevalence': 'Prevalence\n(%)',
              'incidence': 'Incidence\n(per 100 PY)',
              'remission': 'Remission\n(per 100 PY)',
              'excess-mortality': 'Excess Mortality\n(per 1000 PY)' }
    ticks = { 'prevalence': [[.05, .10, .15], [5, 10, 15]],
              'incidence': [[.01, .02, .03], [1, 2, 3]],
              'remission': [[.05, .10, .15], [5, 10, 15]],
              'excess-mortality': [[.001, .01, .1], [1, 10, 100]] }
    ymax = { 'prevalence': .20,
             'incidence': .04,
             'remission': .20,
             'excess-mortality': .20 }
    for ii, type in enumerate(['prevalence', 'incidence', 'remission', 'excess-mortality']):
        # plot ground truth of sim
        data = all_data_hash.get(type, region, year, sex)
        truth = [[d['age_start'], d['truth']] for d in sorted(data, key=lambda dd: dd['age_start'])]

        truth = np.array(truth)

        axes([.3, .15 + (.8*.25)*(3-ii), .6, .8*.25], frameon=False)

        plot(truth[:,0], truth[:,1], color='k',
             linestyle='dotted', linewidth=2, zorder=10)

        # plot noised data for this region
        data = data_hash.get(type, region, year, sex)
        dismod3.plotting.plot_intervals(dm, data, print_sample_size=False,
                                        alpha=.95, color='black',
                                        zorder=900)

        # plot noised data for all regions, sexes, years
        data = data_hash.get(type)
        print type, 'data cnt', len(data)
        dismod3.plotting.plot_intervals(dm, data, print_sample_size=False,
                                        alpha=.75, color='grey', zorder=800)

        # plot fitted value
        dismod3.plotting.plot_mcmc_fit(dm, dismod3.utils.gbd_key_for(type, region, year, sex),
                                       color=dismod3.plotting.color_for[type],
                                       show_data_ui=False)

        # decorate figure
        axis([0,99,0,ymax[type]])
        
        if ii != 3:
            xticks([])
            #xticks([20, 40, 60, 80], ['', '', '', '']) # get gridlines w/o numbers
        else:
            xticks([20, 40, 60, 80])
            xlabel('Age (Years)')
            semilogy([.1])
            axis([0, 99, .0005, .5])

        yticks(*ticks[type])
        figtext(.05, .15 + (.8*.25)*(3-ii+.75), 'abcd'[ii] + ')\n ',
             horizontalalignment='center', verticalalignment='top')
        figtext(.15, .15 + (.8*.25)*(3-ii+.75), label[type], horizontalalignment='center', verticalalignment='top')
        #grid(linestyle='-', color='grey')
        
    savefig('/home/j/Project/dismod/t/sim_%d_%f-posterior.png' % (n, cv))
    savefig('/home/j/Project/dismod/t/sim_%d_%f-posterior.eps' % (n, cv))
    #for k in dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex]):
    #    if dm.vars[k].get('data'):
    #        dismod3.plotting.plot_posterior_predicted_checks(dm, k)
    #        savefig('sim_%d_%f_%s-check.png' % (n, cv, k))

def iqr(X):
    from pymc.utils import hpd
    lb, ub = hpd(X, .4)
    return '(%.0f, %.0f)' % (round(lb-.5), round(ub+.5))

if __name__ == '__main__':
    import optparse
    usage = 'usage: %prog [options] n cv'
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()

    if len(args) == 0:
        # print summary results
        Y = []
        print '\\hline'
        print '\t&\t'.join(['$n$ (rows)', '$cv$ (\\%)', 'repetitions', 'MARE (\\%)','IQR', 'CP (\\%)', 'IQR']) + '\\\\'
        print '\\hline'
        print '\\hline'
        import glob
        for fname in sorted(glob.glob('score*.txt')):
            X = csv2rec(fname, names=['mare', 'coverage'])
            n = float(fname.split('_')[1])
            cv = float(fname.split('_')[2][:-4])

            print '\t&\t'.join(['%d'%n, '%.0f'%cv, '%d'%len(X), '%.1f'%median(X.mare), iqr(X.mare), '%.1f'%median(X.coverage*100), iqr(X.coverage*100)]) + '\\\\'
            Y.append([n, cv, len(X), X.mare.mean(), X.mare.std(), X.coverage.mean()*100, X.coverage.std()*100])
        print '\\hline'
        f = open('/home/j/Project/Models/dm_scores.csv', 'w')
        import csv
        csv_f = csv.writer(f)
        csv_f.writerows([['n', 'cv', 'reps', 'mu_mare', 'std_mare', 'mu_coverage', 'std_coverage']] + Y)
        f.close()
            
    else:
        if len(args) == 2:
            try:
                n = int(args[0])
            except ValueError:
                parser.error('n must be an integer')
            try:
                cv = float(args[1])
                assert cv > 0, 'cv must be positive'
            except ValueError:
                parser.error('cv must be a float')

            dm = fit_simulated_disease(n, cv)

        else:
            n_list = [30,300,3000]
            cv_list = [2, 20, 200]
            for r in range(16):
                for n in n_list:
                    for cv in cv_list:
                        dm = fit_simulated_disease(n, cv)
