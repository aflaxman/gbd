""" Script to distribute the validation jobs on the
IHME cluster"""

import pandas

output_dir = '/home/j/Project/dismod'
validation_name = 'age_group_validation'

def tally_results():
    import glob
    results = pandas.DataFrame(columns='model bias mae mare pc'.split())
    for fname in sorted(glob.glob('%s/%s/*-*.csv' % (output_dir, validation_name))):
        model, rep = fname.split('/')[-1].split('-')

        df = pandas.read_csv(fname, index_col=None)
        df['model'] = model
        results = results.append(df, ignore_index=True)

    summary = pandas.DataFrame()
    for rm, df_rm in results.groupby('model'):
        sum_rm = df_rm.describe()
        sum_rm['index'] = sum_rm.index
        sum_rm['model'] = rm
        summary = summary.append(sum_rm, ignore_index=True)
        
    summary.to_csv('%s/%s/summary_results.csv' % (output_dir, validation_name))
        
    return summary


if __name__ == '__main__':
    results = tally_results()
    print results

