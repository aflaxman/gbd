""" Script to distribute the validation jobs on the
IHME cluster, run each job, and compile the results"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

import optparse
import subprocess

import pandas
import pylab as pl
import pymc as mc
pl.seterr('ignore')
# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import dismod3
reload(dismod3)

import validate_covariates
reload(validate_covariates)

output_dir = '/home/j/Project/dismod'
#output_dir = '/var/tmp/dismod_working'
validation_name = 'fe_model_validation'


def tally_results():
    import glob
    results = pandas.DataFrame(columns='N delta p param bias mae mare pc'.split())
    for fname in sorted(glob.glob('%s/%s/*-*-*-*.csv' % (output_dir, validation_name))):
        N, delta, p, rep = fname.split('/')[-1].split('-')

        df = pandas.read_csv(fname, index_col=None)
        df['N'] = int(N)
        df['delta'] = float(delta)
        df['p'] = float(p)
        results = results.append(df, ignore_index=True)

    results = results.groupby(['N', 'delta', 'p', 'param']).describe()
    results.to_csv('%s/%s/summary_results.csv' % (output_dir, validation_name))
        
    return results

def run_all():
    subprocess.call('mkdir -p %s/%s' % (output_dir, validation_name), shell=True)
    subprocess.call('mkdir -p %s/%s/log/' % (output_dir, validation_name), shell=True)

    names = []
    for N in '10 100 1000 10000'.split():
        for delta in '1 3 9'.split():
            for p in '1 5 10'.split():
                for replicate in range(100):
                    o = '%s/%s/log/%s-%s-%s-%s.txt' % (output_dir, validation_name, N, delta, p, replicate)
                    name_str = '%s-%s-%s-%s-%s' % (validation_name, N, delta, p, replicate)
                    names.append(name_str)

                    call_str = 'qsub -cwd -o %s -e %s ' % (o,o) \
                               + '-N %s ' % name_str \
                               + 'run_on_cluster.sh '

                    call_str += 'tests/run_%s.py --N=%s --delta=%s --p=%s --replicate=%s' % (validation_name, N, delta, p, replicate)

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
    parser.add_option('-N', '--N', default='100',
                      help='rate distribution name')
    parser.add_option('-d', '--delta', default='3',
                      help='rate distribution name')
    parser.add_option('-p', '--p', default='1',
                      help='rate distribution name')
    parser.add_option('-r', '--replicate', default='0',
                      help='replicate number, for saving')
    (options, args) = parser.parse_args()

    if options.runall.lower()=='true':
        run_all()
    elif options.tally.lower()=='true':
        results = tally_results()
        print 'mean over all replicates of median absolute relative error'
        print results
    else:
        N = int(options.N)
        delta = float(options.delta)
        p = int(options.p)
        replicate = int(options.replicate)

        print 'Running validation for:'
        print 'N', options.N
        print 'delta', options.delta
        print 'p', options.p
        print 'replicate', replicate

        model = validate_covariates.validate_covariate_model_fe(N=N, delta_true=delta, beta_true=pl.ones(p), replicate=replicate)
        model.results.to_csv('%s/%s/%s-%s-%s-%s.csv' % (output_dir, validation_name, N, delta, p, replicate))
                             
        

