#!/usr/bin/python2.5
""" Fit continuous single parameter model using gp_re_a model

Example
-------

$ python fit_continuous_spm.py 4773

"""
import matplotlib
matplotlib.use("AGG")

import optparse
import os
import subprocess

import pylab as pl

import dismod3
from dismod3.utils import clean, gbd_keys, type_region_year_sex_from_key


def fit_continuous_spm(id):
    """ Fit continuous single parameter model

    Parameters
    ----------
    id : int
      The model id number for the job to fit

    Example
    -------
    >>> import fit_continuous_spm
    >>> fit_continuous_spm.fit_continuous_spm(4773)
    """

    dm = dismod3.get_disease_model(id)
    
    ## convert model to csv file
    column_names = ['region', 'country', 'year', 'age', 'y', 'se', 'x0', 'x1', 'w0']
    data_list = []

    # add all the model data to the data list
    param_type = 'continuous single parameter'
    for d in dm.filter_data(data_type=param_type):
        row = {}
        row['region'] = dismod3.utils.clean(d['gbd_region'])
        row['country'] = d['country_iso3_code']
        
        row['year'] = int(.5 * (d['year_start'] + d['year_end']))
        row['age'] = int(.5 * (d['age_start'] + d['age_end']))

        row['y'] = d['parameter_value'] * float(d['units'])
        row['se'] = d['standard_error'] * float(d['units'])

        row['x0'] = 1.
        row['x1'] = .1 * (row['year']-1997.)

        row['w0'] = .1 * (row['year']-1997.)

        data_list.append(row)


    # add the time/age/regions that we want to predict to the data list as well
    age_mesh = [0, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    index_dict = {}
    for r in dismod3.gbd_regions[:8]: # FIXME: now i just take a few regions, for fast testing
        for y in [1990, 2005]:
            for a in age_mesh:
                row = {}
                row['region'] = dismod3.utils.clean(r)
                row['country'] = 'all'

                row['year'] = y
                row['age'] = a

                row['y'] = pl.nan
                row['se'] = pl.inf

                row['x0'] = 1.
                row['x1'] = .1 * (row['year']-1997.)

                row['w0'] = .1 * (row['year']-1997.)

                index_dict[(dismod3.utils.clean(r),y,a)] = len(data_list)
                data_list.append(row)


    # save the csv file
    import csv
    fname = dismod3.settings.JOB_WORKING_DIR % id + '/data.csv'

    try:
        f = open(fname, 'w')
        csv.writer(f).writerow(column_names)
        csv.DictWriter(f, column_names).writerows(data_list)
        f.close()
    except IOError, e:
        print 'Warning: could not create data csv.  Maybe it exists already?\n%s' % e


    ## fit the model
    data = pl.csv2rec(fname)
    
    print 'generating model'
    from space_time_model import model
    reload(model)  # for development, automatically reload in case model.py has changed
    mod_mc = model.gp_re_a(data)

    print 'fitting model with mcmc'
    iter = 10000
    #iter = 100 # for testing
    mod_mc.sample(iter, iter/2, 1+iter/2000, verbose=1)
    

    # generate plots of results
    print 'summarizing results'
    param_predicted_stats = mod_mc.param_predicted.stats()
    
    for r in pl.unique(data.region):
        for t in [1990, 2005]:
            x = []
            y = []
            yl = []
            yu = []
            for a in age_mesh:
                x.append(a)
                y.append(param_predicted_stats['mean'][index_dict[(r, t, a)]])
                yl.append(param_predicted_stats['95% HPD interval'][index_dict[(r, t, a)],0])
                yu.append(param_predicted_stats['95% HPD interval'][index_dict[(r, t, a)],1])

            print r, t, zip(x,y)

            key = dismod3.gbd_key_for(param_type, r, t, 'all')
            est = dismod3.utils.interpolate(x, y, dm.get_estimate_age_mesh())
            dm.set_mcmc('mean', key, est)

            est = dismod3.utils.interpolate(x, yl, dm.get_estimate_age_mesh())
            dm.set_mcmc('lower_ui', key, est)

            est = dismod3.utils.interpolate(x, yu, dm.get_estimate_age_mesh())
            dm.set_mcmc('upper_ui', key, est)

            dismod3.tile_plot_disease_model(dm, [key], defaults={})
            try:
                pl.savefig(dismod3.settings.JOB_WORKING_DIR % id + '/dm-%d-posterior-%s-%s-%s.png' % (id, dismod3.utils.clean(r), 'all', t))   # TODO: refactor naming into its own function
            except IOError, e:
                print 'Warning: could not create png.  Maybe it exists already?\n%s' % e

    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    dismod3.try_posting_disease_model(dm, ntries=5)

    print
    print '********************'
    print 'computation complete'
    print '********************'

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

    fit_continuous_spm(id)
        

if __name__ == '__main__':
    main()
