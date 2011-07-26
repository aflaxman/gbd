#!/usr/bin/python2.5
""" Generate map output csv for a given disease model

Example
-------

$ python generate_map_output.py 11094 >map_output.csv

"""

import sys
import csv
import StringIO
import pylab as pl

import dismod3

output_fields = 'region,year,sex,prevalence50,prevalence025,prevalence975,incidence50,incidence025,incidence975,remission50,remission025,remission975,excess-mortality50,excess-mortality025,excess-mortality975,mortality50,mortality025,mortality975,duration50,duration025,duration975'.split(',')

def weighted_average(x, p):
    if len(x) != len(p):
        return None
    return pl.dot(x, p) / pl.sum(p)

def generate_map_output(dm):

    # calculate world population for age-standardization
    age_end = 100
    age_start = 0
    world_pop = pl.zeros(age_end - age_start + 1)
    year = 2005
    sex = 'total'
    for region in dismod3.gbd_regions:
        population_region = dismod3.table.population_by_region_year_sex(dismod3.utils.clean(region), year, sex)[age_start:age_end + 1]
        for age in range(age_end - age_start + 1):
            world_pop[age] += population_region[age]
                                

    # initialize data list
    data = []

    # make first entry of data list the field names
    d = {}
    for f in output_fields:
        d[f] = f
    # fix mortality column heading to read "mortality-with"
    d['mortality50'] = 'mortality-with50'
    d['mortality025'] = 'mortality-with025'
    d['mortality975'] = 'mortality-with975'
    data.append(d)

    # for each region, year, sex append the relevant information
    for r in dismod3.utils.gbd_regions:
        for y in dismod3.utils.gbd_years:
            for s in dismod3.utils.gbd_sexes:
                d = dict(region=r, year=y, sex=s)

                for t in ['prevalence', 'incidence', 'remission', 'excess-mortality', 'mortality', 'duration']:
                    key = dismod3.utils.gbd_key_for(t, r, y, s)

                    for stat, col_name in [['lower_ui', '%s025'%t],
                                           ['upper_ui', '%s975'%t],
                                           ['median', '%s50'%t]]:
                        f = dm.get_mcmc(stat, key)
                        d[col_name] = '%.6f' % weighted_average(f, world_pop)

                data.append(d)
                
    return data

def main():
    argv = sys.argv
    assert len(argv) == 2, 'usage: python generate_map_output.py model_id'

    # download the requested model
    id = int(argv[1])
    dm = dismod3.fetch_disease_model(id)

    # generate a list of dicts of requested data
    map_output_list = generate_map_output(dm)

    # dump list of dicts as csv (in memory, not on disk, hence the StringIO)
    str_io = StringIO.StringIO()
    csv_f = csv.DictWriter(str_io, output_fields)
    csv_f.writerows(map_output_list)

    # print out csv file
    str_io.seek(0)
    print str_io.read()

if __name__ == '__main__':
    main()
