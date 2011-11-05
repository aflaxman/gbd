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

import validate_consistent_model
reload(validate_consistent_model)

output_dir = '/home/j/Project/dismod'
output_dir = '/var/tmp/dismod_working'
validation_name = 'consistent_validation'

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
    results.to_csv('/home/j/Project/dismod/%s/summary_results.csv'%validation_name)
        


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
    parser.add_option('-r', '--replicate', default='1',
                      help='replicate number, for saving')
    (options, args) = parser.parse_args()

    if options.runall.lower()=='true':
        run_all()
    elif options.tally.lower()=='true':
        tally_results()
    else:
        N = int(options.numberofrows)
        delta_true = float(options.delta)
        replicate = int(options.replicate)

        print 'Running random effects validation for:'
        print 'N', N
        print 'delta_true', delta_true
        print 'replicate', replicate

        M = gp.Mean(validate_consistent_model.constant)
        C = gp.Covariance(gp.matern.euclidean, amp=1., diff_degree=2, scale=50)
        gp.observe(M, C, [0, 100], [-5, -5])
        
        true = {}
        for t in 'irf':
            log_rate = gp.Realization(M, C)
            true[t] = lambda x: pl.exp(log_rate(x))

        model = validate_consistent_model.validate_consistent_model_sim(N, delta_true, true)
        model.results.to_csv('%s/%s/%s-%s-%s-%s.csv' % (output_dir, validation_name, options.numberofrows, options.delta, '', options.replicate))
                             
        

