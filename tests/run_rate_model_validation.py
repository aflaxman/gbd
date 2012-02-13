""" Script to distribute the validation jobs on the
IHME cluster, run each job, and compile the results"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

import optparse
import subprocess

import pandas
import dismod3
reload(dismod3)

import validate_rate_model
reload(validate_rate_model)

output_dir = '/home/j/Project/dismod'
#output_dir = '/var/tmp/dismod_working'
validation_name = 'rate_model_validation'

def tally_results():
    import glob
    results = pandas.DataFrame(columns='rate_model bias mae mare pc'.split())
    for fname in sorted(glob.glob('%s/%s/*-*.csv' % (output_dir, validation_name))):
        rate_model, rep = fname.split('/')[-1].split('-')

        df = pandas.read_csv(fname, index_col=None)
        df['rate_model'] = rate_model
        results = results.append(df, ignore_index=True)

    results = results.groupby(['rate_model']).describe()
    results.to_csv('%s/%s/summary_results.csv' % (output_dir, validation_name))
        
    return results

def run_all():
    subprocess.call('mkdir -p %s/%s' % (output_dir, validation_name), shell=True)
    subprocess.call('mkdir -p %s/%s/log/' % (output_dir, validation_name), shell=True)

    names = []
    for rate_type in 'binom beta_binom poisson neg_binom normal log_normal offset_log_normal'.split():
        for replicate in range(1000):
            o = '%s/%s/log/%s-%s.txt' % (output_dir, validation_name, rate_type, replicate)
            name_str = '%s-%s-%s' % (validation_name, rate_type, replicate)
            names.append(name_str)

            call_str = 'qsub -cwd -o %s -e %s ' % (o,o) \
                       + '-N %s ' % name_str \
                       + 'run_on_cluster.sh '

            call_str += 'tests/run_%s.py --ratetype=%s --replicate=%s' % (validation_name, rate_type, replicate)

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
    parser.add_option('-T', '--ratetype', default='neg_binom',
                      help='rate distribution name')
    parser.add_option('-r', '--replicate', default='0',
                      help='replicate number, for saving')
    (options, args) = parser.parse_args()

    if options.runall.lower()=='true':
        run_all()
    elif options.tally.lower()=='true':
        results = tally_results()
        print 'mean over all replicates of median absolute relative error'
        print results.unstack()['mare', 'mean'].unstack()
    else:
        replicate = int(options.replicate)

        print 'Running validation for:'
        print 'ratetype', options.ratetype
        print 'replicate', replicate

        model = validate_rate_model.validate_rate_model(options.ratetype, replicate)
        model.results.to_csv('%s/%s/%s-%s.csv' % (output_dir, validation_name, options.ratetype, options.replicate))
                             
        

