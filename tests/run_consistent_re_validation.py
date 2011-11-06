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

import validate_consistent_re_model
reload(validate_consistent_re_model)

output_dir = '/home/j/Project/dismod'
#output_dir = '/var/tmp/dismod_working'
validation_name = 'consistent_re_validation'

def tally_results():
    import glob
    results = pandas.DataFrame(columns='N delta sigma param bias mae mare pc'.split())
    for fname in sorted(glob.glob('%s/%s/*-*-*-*.csv' % (output_dir, validation_name))):
        N, delta, sigma, rep = fname.split('/')[-1].split('-')

        df = pandas.read_csv(fname, index_col=None)
        df['N'] = int(N)
        df['delta'] = float(delta)
        df['sigma'] = float(sigma)
        results = results.append(df, ignore_index=True)

    results = results.groupby(['N', 'delta', 'sigma', 'param']).describe()
    #results = results.groupby(['N', 'delta', 'sigma', 'param']).mean().drop(['index'], 1)
    results.to_csv('%s/%s/summary_results.csv' % (output_dir, validation_name))
        
    return results

def run_all():
    subprocess.call('mkdir %s/%s' % (output_dir, validation_name), shell=True)
    subprocess.call('mkdir %s/%s/log/' % (output_dir, validation_name), shell=True)

    names = []
    for N in '100 1000 10000'.split():
        for delta in '.1 1'.split():
            for sigma in '.01 .1 1.'.split():
                for replicate in range(12):

                    o = '%s/%s/log/%s-%s-%s-%s.txt' % (output_dir, validation_name, N, delta, sigma, replicate)
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
    parser.add_option('-N', '--numberofrows', default='100',
                      help='number of rows of data to simulate')
    parser.add_option('-d', '--delta', default='.1',
                      help='true over-dispersion parameter delta')
    parser.add_option('-s', '--sigma', default='.5',
                      help='area random effect dispersion sigma_alpha')
    parser.add_option('-r', '--replicate', default='1',
                      help='replicate number, for saving')
    (options, args) = parser.parse_args()

    if options.runall.lower()=='true':
        run_all()
    elif options.tally.lower()=='true':
        results = tally_results()
        print 'mean over all replicates of median absolute relative error'
        print results.unstack()['mare', 'mean'].unstack()
    else:
        N = int(options.numberofrows)
        delta_true = float(options.delta)
        sigma_true = float(options.sigma)*pl.ones(5)
        replicate = int(options.replicate)

        print 'Running random effects validation for:'
        print 'N', N
        print 'delta_true', delta_true
        print 'sigma_true', sigma_true
        print 'replicate', replicate

        M = gp.Mean(validate_consistent_re_model.quadratic)
        C = gp.Covariance(gp.matern.euclidean, amp=1., diff_degree=2, scale=50)
        gp.observe(M, C, [0, 25, 100], [-5, -3, -5])
        
        true = {}
        li = gp.Realization(M, C)
        true['i'] = lambda x: pl.exp(li(x))
        lr = gp.Realization(M, C)
        true['r'] = lambda x: pl.exp(lr(x))
        lf = gp.Realization(M, C)
        true['f'] = lambda x: pl.exp(lf(x))

        model = validate_consistent_re_model.validate_consistent_re(N, delta_true, sigma_true, true)
        model.results.to_csv('%s/%s/%s-%s-%s-%s.csv' % (output_dir, validation_name, options.numberofrows, options.delta, options.sigma, options.replicate))
                             
        

