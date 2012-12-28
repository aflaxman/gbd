""" Script to distribute the validation jobs on the
IHME cluster"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

import optparse
import subprocess

import pandas
import pylab as pl
pl.seterr('ignore')
# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import dismod3
reload(dismod3)

import validate_rate_model
reload(validate_rate_model)

output_dir = '/home/j/Project/dismod'
#output_dir = '/var/tmp/dismod_working'
validation_name = 'rate_model_validation'

def run_all():
    subprocess.call('mkdir -p %s/%s' % (output_dir, validation_name), shell=True)
    subprocess.call('mkdir -p %s/%s/log/' % (output_dir, validation_name), shell=True)

    names = []
    for rate_type in 'poisson neg_binom binom beta_binom normal log_normal offset_log_normal'.split():
        for replicate in range(100):
            o = '%s/%s/log/%s-%s.txt' % (output_dir, validation_name, rate_type, replicate)
            name_str = '%s-%s-%s' % (validation_name, rate_type, replicate)
            names.append(name_str)

            call_str = 'qsub -cwd -o %s -e %s ' % (o,o) \
                       + '-N %s ' % name_str \
                       + 'run_on_cluster.sh '

            call_str += 'tests/validate_rate_model.py %s %d' % (rate_type, replicate)

            print call_str
            subprocess.call(call_str, shell=True)
            
    # after all posteriors have finished running, upload disease model json
    hold_str = '-hold_jid %s ' % ','.join(names)
    o = '%s/%s/log/tally.txt' % (output_dir, validation_name)
    call_str = 'qsub -cwd -o %s -e %s ' % (o,o) \
               + hold_str \
               + '-N %s_tally '%validation_name \
               + 'run_on_cluster.sh '
    call_str += 'tests/tally_%s.py' % validation_name
    subprocess.call(call_str, shell=True)


if __name__ == '__main__':
    run_all()
                             
        

