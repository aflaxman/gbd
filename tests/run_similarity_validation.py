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

pl.seterr('ignore')

def tally_results():
    import glob
    results = pandas.DataFrame(columns='N delta sigma param bias mae mare pc'.split())
    for fname in sorted(glob.glob('%s/%s/*-*-*-*-*-*.csv' % (output_dir, validation_name))):
        N, delta, het, bias, sigma, rep = fname.split('/')[-1].split('-')

        df = pandas.read_csv(fname, index_col=None)
        df['N'] = int(N)
        df['delta'] = float(delta)
        df['heterogeneity'] = het
        df['bias'] = float(bias)
        df['sigma'] = float(sigma)
        results = results.append(df, ignore_index=True)

    results = results.groupby(['N', 'delta', 'heterogeneity', 'bias', 'sigma', 'param']).describe()
    results.to_csv('%s/%s/summary_results.csv' % (output_dir, validation_name))

    return results


def run_all():
    names = []
    bias = .25
    for N in '2 5 10'.split():
        for delta in '.01 1.'.split():
            for sigma in '.125 .5'.split():
                for replicate in range(10):
                    
                    o = '%s/%s/log/%s-%s-%s-%s.txt' % (output_dir, validation_name, N, delta, sigma, replicate)
                    name_str = '%s-%s-%s-%s-%s' % (validation_name, N, delta, sigma, replicate)
                    names.append(name_str)
                
                    call_str = 'qsub -cwd -o %s -e %s ' % (o,o) \
                        + '-N %s ' % name_str \
                        + 'run_on_cluster.sh '

                    #call_str = 'python '
                    call_str += 'tests/run_%s.py -N %s -d %s -r %s -b %s -s %s' % (validation_name, N, delta, replicate, bias, sigma)
                
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
    parser.add_option('-b', '--bias', default='0',
                      help='bias parameter')
    parser.add_option('-s', '--sigma', default='0',
                      help='prior uncertainty sigma parameter')
    parser.add_option('-r', '--replicate', default='0',
                      help='replicate number, for saving')
    (options, args) = parser.parse_args()

    subprocess.call('mkdir %s/%s' % (output_dir, validation_name), shell=True)
    subprocess.call('mkdir %s/%s/log/' % (output_dir, validation_name), shell=True)

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
        replicate = int(options.replicate)
        bias = float(options.bias)
        sigma_prior = float(options.sigma)

        print 'Running random effects validation for:'
        print 'N', N
        print 'delta_true', delta_true
        print 'bias', bias
        print 'sigma_prior', sigma_prior
        print 'replicate', replicate

        M = gp.Mean(validate_similarity.quadratic)
        C = gp.Covariance(gp.matern.euclidean, amp=1., diff_degree=2, scale=50)
        gp.observe(M, C, [0, 30, 100], [-5, -3, -5])

        true = {}
        lp = gp.Realization(M, C)
        true_p = lambda x: pl.exp(lp(x))

        model = validate_similarity.generate_data(N, delta_true, true_p, 'Unusable', bias, sigma_prior)
        for het in 'Very Moderately Slightly'.split():
            model.parameters['p']['heterogeneity'] = het
            validate_similarity.fit(model)
            model.results.to_csv('%s/%s/%s-%s-%s-%s-%s-%s.csv' % (output_dir, validation_name, options.numberofrows, options.delta, het, bias, sigma_prior, options.replicate))
