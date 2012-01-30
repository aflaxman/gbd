#!/usr/bin/python2.5
""" Generate empirical prior of specified parameter type

Expects the disase model json to be saved already.
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

import dismod3
import glob
import pylab as pl
import pymc as mc
import pandas

def upload_fits(id):
    """ Send results of cluster fits to dismod server

    Parameters
    ----------
    id : int
      The model id number

    Example
    -------
    >>> import fit_emp_prior
    >>> fit_emp_prior.fit_emp_prior(2552, 'incidence')
    >>> import upload_fits
    >>> upload_fits.upload_fits(2552)
    """
    # load disease model
    dm = dismod3.load_disease_model(id)  # this merges together results from all fits

    # save dta output
    dir = dismod3.settings.JOB_WORKING_DIR % id  # TODO: refactor into a function
    #dm_to_dta(dm, '%s/regional_predictions' % dir)

    # plot empirical priors (in a separate script, to run after all empirical priors are computed)
    for effect in ['alpha', 'beta', 'gamma', 'delta']:
        try:
            dismod3.plotting.plot_empirical_prior_effects([dm], effect)
            dm.savefig('dm-%d-emp-prior-%s.png' % (id, effect))
        except Exception:
            print 'failed to plot %s' % effect

    # save table output
    try:
        dismod3.table.make_tables(dm)
    except Exception, e:
        print 'Failed to make table'
        print e

    # send to website
    dismod3.try_posting_disease_model(dm, ntries=5)

    # record that job is done
    o = '%s/empirical_priors/stdout/%d_running.txt' % (dir, id)
    f = open(o, 'a')
    import time
    f.write('\n**** JOB DONE AT %s' % time.strftime('%c'))
    f.close()

def merge_data_csvs(id):

    df = pandas.DataFrame()

    dir = dismod3.settings.JOB_WORKING_DIR % id
    #print dir
    for f in sorted(glob.glob('%s/posterior/data-*.csv'%dir)):
        #print 'merging %s' % f
        df2 = pandas.read_csv(f, index_col=None)
        df2.index = df2['index']
        df = df.drop(set(df.index)&set(df2.index)).append(df2)

    df['residual'] = df['value'] - df['mu_pred']
    df['scaled_residual'] = df['residual'] / pl.sqrt(df['value'] * (1 - df['value']) / df['effective_sample_size'])
    #df['scaled_residual'] = df['residual'] * pl.sqrt(df['effective_sample_size'])  # including 
    df['abs_scaled_residual'] = pl.absolute(df['scaled_residual'])

    d = .005 # TODO: save delta in these files, use negative binomial to calc logp
    df['logp'] = [mc.negative_binomial_like(x*n, (p+1e-3)*n, d*(p+1e-3)*n) for x,p,n in zip(df['value'], df['mu_pred'], df['effective_sample_size'])]
    df['logp'][df['data_type'] == 'rr'] = df['scaled_residual'][df['data_type'] == 'rr']

    df = df.sort('logp')

    #print df.filter('data_type area age_start age_end year_start sex effective_sample_size value residual logp'.split())[:25]
    return df

import csv, subprocess
population_by_age = dict(
    [[(dismod3.utils.clean(r['Country Code']), int(r['Year']), r['Sex']),
      [max(.001,float(r['Age %d Population' % i])) for i in range(dismod3.settings.MAX_AGE)]] 
     for r in csv.DictReader(open(dismod3.settings.CSV_PATH + 'population.csv'))
     ]
)
def dm_to_dta(dm, fname):
    X = ['type, region, sex, year, age, pop, prior, posterior, upper, lower'.split(', ')]

    for t in dismod3.utils.output_data_types:
        for r in dismod3.settings.gbd_regions:
            r = dismod3.utils.clean(r)
            for s in ['male', 'female']:
                for y in [1990, 2005, 2010]:
                    k = dismod3.utils.gbd_key_for(t, r, y, s)

                    prior = dm.get_mcmc('emp_prior_mean', k)
                    if len(prior) == 0:
                        prior = -99 * pl.ones(100)

                    posterior = dm.get_mcmc('mean', k)
                    lower = dm.get_mcmc('lower_ui', k)
                    upper = dm.get_mcmc('upper_ui', k)
                    if len(posterior) == 0:
                        posterior = -99 * pl.ones(100)
                        lower = -99 * pl.ones(100)
                        upper = -99 * pl.ones(100)
                    for a in range(100):
                        X.append([t, r, s, y, a,
                                  population_by_age[r,y,s][a],
                                  prior[a],
                                  posterior[a],
                                  upper[a],
                                  lower[a]
                                  ])

            f = open('%s.csv'%fname, 'w')
            csv.writer(f).writerows(X)
            f.close()

            convert_cmd = 'echo \'library(foreign); X=read.csv("%s"); write.dta(X, "%s")\' | %s --no-save' % ('%s.csv'%fname, '%s.dta'%fname, dismod3.settings.R_PATH)
    ret = subprocess.call(convert_cmd, shell=True)
    assert ret == 0, 'return code %d' % ret


def main():
    import optparse

    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    upload_fits(id)
      

if __name__ == '__main__':
    main()
