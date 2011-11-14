""" Script to distribute the validation jobs on the
IHME cluster, run each job, and compile the results"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

import optparse
import subprocess

import pylab as pl
import pandas
from pymc import gp

import validate_age_integrating_re
reload(validate_age_integrating_re)

output_dir = '/home/j/Project/dismod'
#output_dir = '/var/tmp/dismod_working'
validation_name = 'ai_re_validation'

def tally_results():
    import glob
    results = pandas.DataFrame(columns='N delta sigma param bias mae mare pc'.split())
    for fname in sorted(glob.glob('/home/j/Project/dismod/%s/*-*-*-*.csv'%validation_name)):
        N, delta, sigma, rep = fname.split('/')[-1].split('-')

        df = pandas.read_csv(fname, index_col=None)
        df['N'] = int(N)
        df['delta'] = float(delta)
        df['sigma'] = float(sigma)
        results = results.append(df, ignore_index=True)

    results = results.groupby(['N', 'delta', 'sigma', 'param']).describe()
    #results = results.groupby(['N', 'delta', 'sigma', 'param']).mean().drop(['index'], 1)
    results.to_csv('/home/j/Project/dismod/%s/summary_results.csv'%validation_name)
        
    return results

def run_all():
    names = []
    for N in '10 100 1000 10000'.split():
        for delta in '.01 .1 1.'.split():
            for sigma in '.1 .5 2.5'.split():
                for replicate in range(10):

                    o = '/home/j/Project/dismod/%s/log/%s-%s-%s-%s.txt' % (validation_name, N, delta, sigma, replicate)
                    name_str = '%s-%s-%s-%s-%s' % (validation_name, N, delta, sigma, replicate)
                    names.append(name_str)

                    call_str = 'qsub -cwd -o %s -e %s ' % (o,o) \
                        + '-N %s ' % name_str \
                        + 'run_on_cluster.sh '

                    call_str += 'tests/run_%s.py -N %s -d %s -s %s -r %s' % (validation_name, N, delta, sigma, replicate)

                    print call_str
                    subprocess.call(call_str, shell=True)
                    
    # after all posteriors have finished running, upload disease model json
    hold_str = '-hold_jid %s ' % ','.join(names)
    o = '/home/j/Project/dismod/%s/log/tally.txt' % validation_name
    call_str = 'qsub -cwd -o %s -e %s ' % (o,o) \
               + hold_str \
               + '-N %s_tally '%validation_name \
               + 'run_on_cluster.sh '
    call_str += 'tests/run_%s.py --tally=true' % validation_name
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
    parser.add_option('-S', '--smoothing', default='Slightly',
                      help='replicate number, for saving')
    parser.add_option('-r', '--replicate', default='1',
                      help='replicate number, for saving')
    (options, args) = parser.parse_args()

    if options.runall.lower()=='true':
        run_all()
    elif options.tally.lower()=='true':
        results = tally_results()
        print 'count of all replicates'
        print results.unstack()['mare', 'count'].unstack()
        print "\ninspect with:\nresults.unstack()['mare', '50%'].unstack() # for example"
        print "or: results.unstack()['mare', '50%'].unstack(2).reindex(columns='Very Moderately Slightly'.split())"
                                
    else:
        N = int(options.numberofrows)
        delta_true = float(options.delta)
        sigma_true = float(options.sigma)*pl.ones(5)
        replicate = int(options.replicate)
        smoothness = options.smoothing

        print 'Running random effects validation for:'
        print 'N', N
        print 'delta_true', delta_true
        print 'sigma_true', sigma_true
        print 'replicate', replicate
        print 'smoothness', smoothness

        M = gp.Mean(validate_age_integrating_re.quadratic)
        C = gp.Covariance(gp.matern.euclidean, amp=1., diff_degree=2, scale=50)
        gp.observe(M, C, [0, 25, 100], [-5, -3, -5])
        
        log_p = gp.Realization(M, C)
        true_p = lambda x: pl.exp(log_p(x))

        model = validate_age_integrating_re.validate_ai_re(N, delta_true, sigma_true, true_p, smoothness)
        model.results.to_csv('%s/%s/%s-%s-%s-%s.csv' % (output_dir, validation_name, options.numberofrows, options.delta, options.sigma, options.replicate))
                             
        

