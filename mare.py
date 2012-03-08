""" Find Median Absolute Relative Error of a DisMod model
Example
-------

$ ./run_on_cluster.sh mare.py 27944

"""

import optparse
import pylab as pl

import upload_fits
reload(upload_fits)

def mare(id):
    df = upload_fits.merge_data_csvs(id)
    df['are'] = pl.absolute((df['value'] - df['mu_pred']) / df['value'])

    print 'mare by type:'
    print pl.sort(df.groupby('data_type')['are'].median())
    
    print
    print 'overall mare: %.3f' % df['are'].median()

    return df['are'].median()

if __name__ == '__main__':
    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    mare(id)

