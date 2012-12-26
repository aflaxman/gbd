""" Script to distribute the Random Effects validation jobs on the
IHME cluster, run each job, and compile the results"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

import optparse
import subprocess

import pylab as pl
import pandas

pl.seterr('ignore')

import validate_covariates
reload(validate_covariates)


def tally_results():
    import glob
    results = pandas.DataFrame(columns='N delta sigma param bias mae mare pc'.split())
    for fname in sorted(glob.glob('/home/j/Project/dismod/re_validation/*-*-*-*-*.csv')):
        N, delta, sigma, ess, rep = fname.split('/')[-1].split('-')

        df = pandas.read_csv(fname, index_col=None)
        df['N'] = int(N)
        df['delta'] = float(delta)
        df['sigma'] = float(sigma)
        df['ess'] = float(ess)
        results = results.append(df, ignore_index=True)

    results = results.groupby(['N', 'delta', 'sigma', 'ess', 'param']).describe()
    results.to_csv('/home/j/Project/dismod/re_validation/summary_results.csv')
    return results


def run_all():
    names = []
    for N in '10 100'.split():
        for delta in '.001 .01 .1'.split():
            for sigma in '.5'.split():
                for ess in '10 1000'.split():
                    for replicate in range(10):

                        o = '/home/j/Project/dismod/re_validation/log/%s-%s-%s-%s-%s.txt' % (N, delta, sigma, ess, replicate)
                        name_str = 're%s-%s-%s-%s-%s' % (N, delta, sigma, ess, replicate)
                        names.append(name_str)

                        call_str = 'qsub -cwd -o %s -e %s ' % (o,o) \
                            + '-N %s ' % name_str \
                            + 'run_on_cluster.sh '

                        call_str += 'tests/run_re_validation.py -N %s -d %s -s %s -e %s -r %s' % (N, delta, sigma, ess, replicate)

                        print call_str
                        subprocess.call(call_str, shell=True)
                    
    # after all posteriors have finished running, upload disease model json
    hold_str = '-hold_jid %s ' % ','.join(names)
    o = '/home/j/Project/dismod/re_validation/log/tally.txt'
    call_str = 'qsub -cwd -o %s -e %s ' % (o,o) \
               + hold_str \
               + '-N re_tally ' \
               + 'run_on_cluster.sh '
    call_str += 'tests/run_re_validation --tally=true'
    subprocess.call(call_str, shell=True)


if __name__ == '__main__':
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-R', '--runall', default='false',
                      help='run all models on the cluster')
    parser.add_option('-t', '--tally', default='false',
                      help='tally the results')
    parser.add_option('-N', '--numberofrows', default='10',
                      help='number of rows of data to simulate')
    parser.add_option('-d', '--delta', default='.01',
                      help='true over-dispersion parameter delta')
    parser.add_option('-s', '--sigma', default='.1',
                      help='area random effect dispersion sigma_alpha')
    parser.add_option('-e', '--ess', default='100',
                      help='effective sample size')
    parser.add_option('-r', '--replicate', default='1',
                      help='replicate number, for saving')
    (options, args) = parser.parse_args()

    if options.runall.lower()=='true':
        run_all()
    elif options.tally.lower()=='true':
        results = tally_results()
        print 'count of replicates by parameter settings'
        print results.unstack()['mare', 'count'].unstack()
        print 
    else:
        N = int(options.numberofrows)
        delta_true = float(options.delta)
        sigma_true = float(options.sigma)*pl.ones(5)
        ess = float(options.ess)
        replicate = int(options.replicate)

        print 'Running random effects validation for:'
        print 'N', N
        print 'delta_true', delta_true
        print 'sigma_true', sigma_true
        print 'replicate', replicate
        
        model = validate_covariates.validate_covariate_model_re(N, delta_true=delta_true, sigma_true=sigma_true, ess=ess)
        model.results.to_csv('/home/j/Project/dismod/re_validation/%s-%s-%s-%s.csv' % (options.numberofrows, options.delta, options.sigma, options.replicate))
                             
        

