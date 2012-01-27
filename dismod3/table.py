""" All methods for reporting Disease Models in tables
"""

from xlwt import *
import pylab as pl

import dismod3
from dismod3.neg_binom_model import countries_for, population_by_age

def make_tables(dm):
    """Make a table representation of the disease model data and
    estimates provided for a region, a year, a sex

    Parameters
    ----------
    dm : DiseaseJson object
    """
    print 'Making tables'
    wb = Workbook()
    make_count_table(dm, wb)
    make_hazard_table(dm, wb)
    make_population_rate_table(dm, wb)

    make_data_page(dm, wb)
    make_priors_and_covariates_page(dm, wb)

    dir = dismod3.settings.JOB_WORKING_DIR % dm.id
    fname = 'dm-%d.xls' % dm.id
    wb.save('%s/%s'  % (dir, fname))

def make_count_table(dm, wb):
    """Make a table representation of disease counts

    Parameters
    ----------
    dm : DiseaseJson object
    wb : the xlwt Workbook to add the count table to
    """

    ws = wb.add_sheet('Counts')

    params = 'incidence prevalence excess-mortality incidence_x_duration'.split()
    cols = 'region year sex age param mean lower upper'.split()
    for col, col_head in enumerate(cols):
        ws.write(0, col, col_head)

    row = 1
    for region in dismod3.settings.gbd_regions:
        for year in dismod3.settings.gbd_years:
            for sex in dismod3.settings.gbd_sexes:
                for param in params:

                    mean_array = dm.get_mcmc('mean', dismod3.utils.gbd_key_for(param, region, year, sex))
                    lower_array = dm.get_mcmc('lower_ui', dismod3.utils.gbd_key_for(param, region, year, sex))
                    upper_array = dm.get_mcmc('upper_ui', dismod3.utils.gbd_key_for(param, region, year, sex))

                    if len(mean_array) == 0:
                        continue

                    pop = population_by_region_year_sex(region, year, sex)
                    if param == 'excess-mortality':
                        # use estimate of prevalent population
                        pop *= dm.get_mcmc('mean', dismod3.utils.gbd_key_for('prevalence', region, year, sex))
                    elif param == 'incidence_x_duration':
                        pop = pl.ones_like(pop)  # population already included in iX stoch

                    for a in range(len(dismod3.settings.gbd_ages)-1):
                        age = dismod3.settings.gbd_ages[a]

                        ages = range(age, dismod3.settings.gbd_ages[a+1])
                        pop_weights = pop[ages]

                        mean = dismod3.utils.rate_for_range(mean_array, ages, pop_weights)
                        lower = dismod3.utils.rate_for_range(lower_array, ages, pop_weights)
                        upper = dismod3.utils.rate_for_range(upper_array, ages, pop_weights)

                        for col, val in enumerate(cols):
                            ws.write(row, col, vars()[val])

                        row += 1

def make_hazard_table(dm, wb):
    """Make a table representation of disease hazards

    Parameters
    ----------
    dm : DiseaseJson object
    wb : the xlwt Workbook to add the count table to
    """

    ws = wb.add_sheet('Hazards')

    params = 'prevalence incidence remission excess-mortality duration m_with relative-risk'.split()
    cols = 'region year sex age param mean lower upper'.split()
    for col, col_head in enumerate(cols):
        ws.write(0, col, col_head)

    row = 1
    for region in dismod3.settings.gbd_regions:
        for year in dismod3.settings.gbd_years:
            for sex in dismod3.settings.gbd_sexes:
                for param in params:

                    mean_array = dm.get_mcmc('mean', dismod3.utils.gbd_key_for(param, region, year, sex))
                    lower_array = dm.get_mcmc('lower_ui', dismod3.utils.gbd_key_for(param, region, year, sex))
                    upper_array = dm.get_mcmc('upper_ui', dismod3.utils.gbd_key_for(param, region, year, sex))

                    if len(mean_array) == 0:
                        continue

                    pop = population_by_region_year_sex(region, year, sex)
                    if param in ['remission', 'excess-mortality', 'duration', 'm_with']:
                        # use estimate of prevalent population
                        pop *= dm.get_mcmc('mean', dismod3.utils.gbd_key_for('prevalence', region, year, sex))

                    for a in range(len(dismod3.settings.gbd_ages)-1):
                        age = dismod3.settings.gbd_ages[a]

                        ages = range(age, dismod3.settings.gbd_ages[a+1])
                        pop_weights = pop[ages]/pop[ages].sum()  # make weights sum to one

                        mean = dismod3.utils.rate_for_range(mean_array, ages, pop_weights)
                        lower = dismod3.utils.rate_for_range(lower_array, ages, pop_weights)
                        upper = dismod3.utils.rate_for_range(upper_array, ages, pop_weights)

                        for col, val in enumerate(cols):
                            ws.write(row, col, vars()[val])

                        row += 1


def make_population_rate_table(dm, wb):
    """Make a table representation of disease population rates

    Parameters
    ----------
    dm : DiseaseJson object
    wb : the xlwt Workbook to add the count table to
    """

    ws = wb.add_sheet('Population Rates')

    params = 'incidence excess-mortality m_with'.split()
    cols = 'region year sex age param mean lower upper'.split()
    for col, col_head in enumerate(cols):
        ws.write(0, col, col_head)

    row = 1
    for region in dismod3.settings.gbd_regions:
        for year in dismod3.settings.gbd_years:
            for sex in dismod3.settings.gbd_sexes:
                for param in params:

                    mean_array = dm.get_mcmc('mean', dismod3.utils.gbd_key_for(param, region, year, sex))
                    lower_array = dm.get_mcmc('lower_ui', dismod3.utils.gbd_key_for(param, region, year, sex))
                    upper_array = dm.get_mcmc('upper_ui', dismod3.utils.gbd_key_for(param, region, year, sex))

                    if len(mean_array) == 0:
                        continue

                    pop = population_by_region_year_sex(region, year, sex)
                    prev_array = pl.array(dm.get_mcmc('mean', dismod3.utils.gbd_key_for('prevalence', region, year, sex)))

                    if param == 'incidence':
                        mean_array *= 1-prev_array
                        lower_array *= 1-prev_array
                        upper_array *= 1-prev_array
                    else:
                        mean_array *= prev_array
                        lower_array *= prev_array
                        upper_array *= prev_array

                    for a in range(len(dismod3.settings.gbd_ages)-1):
                        age = dismod3.settings.gbd_ages[a]

                        ages = range(age, dismod3.settings.gbd_ages[a+1])
                        pop_weights = pop[ages]/pop[ages].sum()  # make weights sum to one
                        mean = dismod3.utils.rate_for_range(mean_array, ages, pop_weights)
                        lower = dismod3.utils.rate_for_range(lower_array, ages, pop_weights)
                        upper = dismod3.utils.rate_for_range(upper_array, ages, pop_weights)

                        for col, val in enumerate(cols):
                            ws.write(row, col, vars()[val])

                        row += 1


def make_data_page(dm, wb):
    """ Write data as a table that can be loaded into dismod"""
     # don't include all-cause mortality data in data table
    data_list = [d for d in dm.data \
                 if d['data_type'] != 'all-cause mortality data']

    
    ws = wb.add_sheet('data')

    if len(data_list) == 0:
        return

    all_keys = set()
    for d in data_list:
        all_keys |= set(d.keys())

    required_keys = ['GBD Cause', 'Parameter', 'GBD Region', 'Country ISO3 Code',
                     'Sex', 'Year Start', 'Year End', 'Age Start', 'Age End',
                     'Parameter Value', 'Standard Error', 'Units', ]

    redundant_keys = ['_row', 'age_weights', 'id', 'value', 'condition', 'data_type', 'region']

    additional_keys = sorted(all_keys - set([dismod3.utils.clean(k) for k in required_keys] + redundant_keys))

    keys = required_keys + additional_keys

    if len(keys) > 256:  # limitation in old excel format
        ws.write(0, 0, 'could not write data: too many columns')
        return

    
    for c, k in enumerate(keys):
        if k == 'GBD Region':
            k = 'Region'
        ws.write(0, c, k)
    for r, d in enumerate(sorted(data_list, key=lambda d: d.get('_row'))):
        for c, k in enumerate(keys):
            val = d.get(dismod3.utils.clean(k), '')
            if val == 'mortality data':
                val = 'with condition mortality data'
            ws.write(r+1, c, val)
            
def make_priors_and_covariates_page(dm, wb):
    """ Write json for the priors and covariates in the workbook, to ensure the results reproducible"""

    ws = wb.add_sheet('priors')
    for r, type in enumerate(['prevalence', 'incidence', 'remission', 'excess-mortality', 'duration', 'relative-risk']):
        prior_str = dm.get_global_priors(type)
        global_priors = dm.get_global_priors(type)
        ws.write(r, 0, type)
        ws.write(r, 1, prior_str)

    r += 2
    ws.write(r, 0, 'covariate')
    ws.write(r, 1, 'reference value')
    covariates_dict = dm.get_covariates()
    for level in ['Study_level', 'Country_level']:
        for k in covariates_dict.get(level, []):
            if k != 'none' and covariates_dict[level][k]['rate']['value']:
                ref_val = covariates_dict[level][k]['value']['value']

                r += 1
                ws.write(r, 0, k)
                ws.write(r, 1, ref_val)
                

def population_by_region_year_sex(region, year, sex):
    """ Calculate population of specified region, year and sex

    Parameters:
    -----------
    region : str
    year : str
    sex : str

    Return:
    -------
    population : list of float numbers
    """
    region = dismod3.utils.clean(region)
    year = str(year)
    sex = dismod3.utils.clean(sex)
    population = pl.zeros(dismod3.settings.MAX_AGE)
    for iso in countries_for[region]:
        if sex == 'all' or sex == 'total':
            population += population_by_age[(iso, year, 'male')]
            population += population_by_age[(iso, year, 'female')]
        else:
            population += population_by_age[(iso, year, sex)]
    return population

