#!/usr/bin/python2.5
""" Fit disease model using DisMod 4!

Example
-------

$ python fit_dm4.py 4773

"""
import matplotlib
matplotlib.use("AGG")

import optparse
import os
import subprocess

import pylab as pl

import dismod3
from dismod3.utils import clean, gbd_keys, type_region_year_sex_from_key

inf = 1e19

def fit(id):
    """ Download model, conduct fit, and upload results

    Parameters
    ----------
    id : int
      The model id number for the job to fit

Commandline Version:

[omak] dismod4.abie ] test/parameter.sh
[omak] dismod4.abie ] example/simulate.py 5 1 100
[omak] dismod4.abie ] /tmp/dismod4_csv test/parameter.csv measure_in.csv sfun_in.csv
dismod4_csv: Attempt to overwrite the existing file
sfun_in.csv
[omak] dismod4.abie ] /tmp/dismod4_csv test/parameter.csv measure_in.csv sfun_in.csv sfun_out.csv measure_out.csv

    """

    dm = dismod3.get_disease_model(id)
    mort = dismod3.fetch_disease_model('all-cause_mortality')
    dm.data += mort.data

    ## convert model to csv file
    column_names = 'time_lower,time_upper,age_lower,age_upper,likelihood_name,likelihood_sigma,likelihood_beta,value,integrand'.split(',')
    data_list = []

    # add all the model data to the data list
    for d in dm.data:
        row = {}
        row['time_lower'] = d['year_start']
        row['time_upper'] = d['year_end']  # TODO: determine if this should be +1

        row['age_lower'] = d['age_start']+1.
        row['age_upper'] = d['age_end']+1.  # TODO: determine if this should be +1


        row['likelihood_name'] = 'gaussian'
        row['likelihood_sigma'] = .0001  # TODO: use more accurate sigma
        row['likelihood_beta'] = 1.

        row['value'] = d['value'] / float(d.get('units', 1.))

        for dm3_type, dm4_type in [['remission data', 'remission'],
                                   ['excess-mortality data', 'excess'],
                                   ['incidence data', 'incidence'],
                                   ['mrr data', 'risk'],
                                   ['prevalence data', 'prevalence'],
                                   ['all-cause mortality data', 'all_cause'],
                                   ]:
            if d['data_type'] == dm3_type:
                row['integrand'] = dm4_type
                data_list.append(row)
                break

    # add the time/age/regions that we want to predict to the data list as well
    age_mesh = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    index_dict = {}
    for year in [1990, 2005]:
        for age in age_mesh:
            for type in ['remission', 'excess', 'incidence', 'risk', 'prevalence']:
                row = {}

                row['time_lower'] = year
                row['time_upper'] = year

                row['age_lower'] = age
                row['age_upper'] = age

                row['likelihood_name'] = 'gaussian'
                row['likelihood_sigma'] = inf
                row['likelihood_beta'] = 1.

                row['value'] = 0.

                row['integrand'] = type
                
                index_dict[(type, year, age)] = len(data_list)
                data_list.append(row)


    # save the csv file
    import csv
    fname = dismod3.settings.JOB_WORKING_DIR % id + '/measure_in.csv'

    try:
        f = open(fname, 'w')
        csv.writer(f).writerow(column_names)
        csv.DictWriter(f, column_names).writerows(data_list)
        f.close()
    except IOError, e:
        print 'Warning: could not create data csv.  Maybe it exists already?\n%s' % e


    ## fit the model
    dir = dismod3.settings.JOB_WORKING_DIR % id
    # rm %s/sfun_in.csv # if you want to regenerate default mesh parameters
    call_str = '/tmp/dismod4.abie/build/src/dismod4_csv /tmp/dismod4.abie/test/parameter.csv %s/measure_in.csv %s/sfun_in.csv' % (dir, dir)
    subprocess.call(call_str, shell=True)

    call_str = '/tmp/dismod4.abie/build/src/dismod4_csv /tmp/dismod4.abie/test/parameter.csv %s/measure_in.csv %s/sfun_in.csv %s/sfun_out.csv %s/measure_out.csv' % (dir, dir, dir, dir)
    subprocess.call(call_str, shell=True)
    

    # generate plots of results
    print 'summarizing results'
    measure_out = pl.csv2rec('%s/measure_out.csv' % dir)

    r = 'asia_southeast'
    for year in [1990]:
        for sex in ['male']:
            for dm3_type, dm4_type in [['remission', 'remission'],
                                       ['excess-mortality', 'excess'],
                                       ['incidence', 'incidence'],
                                       ['mrr', 'risk'],
                                       ['prevalence', 'prevalence'],
                                       ]:
                x = [0]
                y = [0]
                for age in age_mesh:
                    x.append(age)
                    y.append(measure_out.model[index_dict[(dm4_type, year, age)]])

                key = dismod3.gbd_key_for(dm3_type, r, year, sex)
                est = dismod3.utils.interpolate(x, y, dm.get_estimate_age_mesh())
                dm.set_truth(key, est)

                dismod3.tile_plot_disease_model(dm, [key], defaults={})
                try:
                    pl.savefig(dismod3.settings.JOB_WORKING_DIR % id + '/dm-%d-posterior-%s-%s-%s.png' % (id, dm3_type, sex, year))   # TODO: refactor naming into its own function
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

    fit(id)
        

if __name__ == '__main__':
    main()
