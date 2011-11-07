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

import validate_similarity
reload(validate_similarity)

output_dir = '/home/j/Project/dismod'
#output_dir = '/var/tmp/dismod_working'
validation_name = 'similarity_validation'

def tally_results(output_dir, validation_name):
    import glob
    results = pandas.DataFrame(columns='N delta sigma param bias mae mare pc'.split())
    for fname in sorted(glob.glob('%s/%s/*-*-*-*.csv' % (output_dir, validation_name))):
        N, delta, het, rep = fname.split('/')[-1].split('-')

        df = pandas.read_csv(fname, index_col=None)
        df['N'] = int(N)
        df['delta'] = float(delta)
        results = results.append(df, ignore_index=True)

    results = results.groupby(['N', 'delta', 'heterogeneity', 'param']).describe()
    results.to_csv('%s/%s/summary_results.csv' % (output_dir, validation_name))

    return results


def run_all(output_dir, validation_name):
    names = []
    for N in '1 5 10 50'.split():
        for delta in '.001 .01 .1'.split():
            for het in 'Very Moderately Slightly'.split():
                for replicate in range(10):
                    
                    o = '%s/%s/log/%s-%s-%s-%s.txt' % (output_dir, validation_name, N, delta, het, replicate)
                    name_str = '%s-%s-%s-%s-%s' % (validation_name, N, delta, het, replicate)
                    names.append(name_str)
                
                    call_str = 'qsub -cwd -o %s -e %s ' % (o,o) \
                        + '-N %s ' % name_str \
                        + 'run_on_cluster.sh '
                
                    call_str += 'tests/run_%s.py -N %s -d %s -r %s -H %s' % (validation_name, N, delta, replicate, het)
                
                    print call_str
                    subprocess.call(call_str, shell=True)
                    
    # after all posteriors have finished running, upload disease model json
    hold_str = '-hold_jid %s ' % ','.join(names)
    o = '%s/%s/log/tally.txt' % (output_dir, validation_name)
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
    parser.add_option('-d', '--delta', default='.1',
                      help='true over-dispersion parameter delta')
    parser.add_option('-H', '--heterogeneity', default='Very',
                      help='expert prior on over-dispersion parameter delta (Very, Moderately, Slightly)')
    parser.add_option('-r', '--replicate', default='0',
                      help='replicate number, for saving')
    (options, args) = parser.parse_args()

    subprocess.call('mkdir %s/%s' % (output_dir, validation_name), shell=True)
    subprocess.call('mkdir %s/%s/log/' % (output_dir, validation_name), shell=True)

    if options.runall.lower()=='true':
        run_all()
    elif options.tally.lower()=='true':
        results = tally_results()
        print 'median over all replicates of median absolute relative error'
        print results.unstack()['mare', '50%'].unstack()
        
    else:
        N = int(options.numberofrows)
        delta_true = float(options.delta)
        replicate = int(options.replicate)

        print 'Running random effects validation for:'
        print 'N', N
        print 'delta_true', delta_true
        print 'delta_expert', options.heterogeneity
        print 'replicate', replicate

        M = gp.Mean(validate_similarity.quadratic)
        C = gp.Covariance(gp.matern.euclidean, amp=1., diff_degree=2, scale=50)
        gp.observe(M, C, [0, 30, 100], [-5, -3, -5])

        true = {}
        lp = gp.Realization(M, C)
        true_p = lambda x: pl.exp(lp(x))

        model = validate_similarity.validate_similarity(N, delta_true, true_p, options.heterogeneity)
        model.results.to_csv('%s/%s/%s-%s-%s-%s.csv' % (output_dir, validation_name, options.numberofrows, options.delta, options.heterogeneity, options.replicate))
