""" Script to compile the results of validation jobs run on the
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

def tally_results():
    import glob
    results = pandas.DataFrame(columns='rate_model bias mae mare pc'.split())
    for fname in sorted(glob.glob('%s/%s/*-*-*.csv' % (output_dir, validation_name))):
        rate_model, data_model, rep = fname.split('/')[-1].split('-')

        df = pandas.read_csv(fname, index_col=None)
        df['rate_model'] = rate_model
        df['data_model'] = data_model
        results = results.append(df, ignore_index=True)

    summary = results.groupby(['rate_model', 'data_model']).describe()
    summary.to_csv('%s/%s/summary_results.csv' % (output_dir, validation_name))
        
    return summary

if __name__ == '__main__':
    results = tally_results()
    print results

