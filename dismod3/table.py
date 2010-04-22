""" All methods for reporting Disease Models in tables

Useful High-level Methods::

    table_by_region_year_sex(dm_json, keys, user, group_size)
    table_by_region(dm_json, keys, user, group_size)
    table(dm_json, keys, user, group_size)

Useful Low-level Methods::

    table_region_sheet(dm, keys, wb, name, user, group_size)
    table_disease_model(dm, keys, ws, x, y, group_size)
    write_table_group_value(dm, key, item, ws, x, y, group_sizes)
"""

from xlwt import *
from StringIO import StringIO
from time import strftime
from dismod3.plotting import GBDDataHash 
from dismod3.neg_binom_model import countries_for, population_by_age
import dismod3
from dismod3.utils import clean, rate_for_range
from disease_json import *

def table_by_region_year_sex(dm_json, keys, user, group_size):
    """Make a table representation of the disease model data and
    estimates provided for a region, a year, a sex

    Parameters
    ----------
    dm_json : str or DiseaseJson object
      the json string or a thin python wrapper around this data that
      is to be plotted
    keys : list the keys to include
    user : username
    group_size : positive integer smaller than 102
    """
    if isinstance(dm_json, DiseaseJson):
        dm = dm_json
    else:
        try:
            dm = DiseaseJson(dm_json)
        except ValueError:
            print 'ERROR: dm_json is not a DiseaseJson object or json string'
            return

    keys = [k for k in keys if k.split(KEY_DELIM_CHAR)[0] != 'bins']
    type, region, year, sex = keys[0].split(dismod3.utils.KEY_DELIM_CHAR)
    wb = Workbook()
    ws = wb.add_sheet("%s-%s-%s" % (region, year, sex))
    date = strftime("%Y/%m/%d")
    time = strftime("%H:%M:%S")
    ws.write(0, 0, "Dismod III output, date: %s, time: %s, user: %s" % (date, time, user))
    table_disease_model(dm, keys, ws, 0, 0, group_size)
    f = StringIO()
    wb.save(f)
    f.seek(0)
    return f.read()

def table_by_region(dm_json, keys, user, group_size):
    """Make a table representation of the disease model data and
    estimates provided for a region

    Parameters
    ----------
    dm_json : str or DiseaseJson object
      the json string or a thin python wrapper around this data that
      is to be plotted
    keys : list the keys to include
    user : username
    group_size : positive integer smaller than 102
    """
    if isinstance(dm_json, DiseaseJson):
        dm = dm_json
    else:
        try:
            dm = DiseaseJson(dm_json)
        except ValueError:
            print 'ERROR: dm_json is not a DiseaseJson object or json string'
            return

    keys = [k for k in keys if k.split(KEY_DELIM_CHAR)[0] != 'bins']
    keys_male_1 = []
    keys_male_2 = []
    keys_female_1 = []
    keys_female_2 = []
    y = keys[0].split(KEY_DELIM_CHAR)[2]
    for k in keys:
        if k.split(KEY_DELIM_CHAR)[3] == 'male':
            if k.split(KEY_DELIM_CHAR)[2]  == y:
                keys_male_1.append(k)
            else:
                keys_male_2.append(k)
        else:
            if k.split(KEY_DELIM_CHAR)[2]  == y:
                keys_female_1.append(k)
            else:
                keys_female_2.append(k)

    wb = Workbook()
    table_region_sheet(dm, keys, wb, keys[0].split(KEY_DELIM_CHAR)[1], user, group_size)
    f = StringIO()
    wb.save(f)
    f.seek(0)
    return f.read()

def table(dm_json, keys, user, group_size):
    """Make a work book in which each sheet for a region

    Parameters
    ----------
    dm_json : str or DiseaseJson object
      the json string or a thin python wrapper around this data that
      is to be plotted
    keys : list the keys to include
    user : username
    group_size : positive integer smaller than 102
    """
    if isinstance(dm_json, DiseaseJson):
        dm = dm_json
    else:
        try:
            dm = DiseaseJson(dm_json)
        except ValueError:
            print 'ERROR: dm_json is not a DiseaseJson object or json string'
            return

    keys = [k for k in keys if k.split(KEY_DELIM_CHAR)[0] != 'bins'] 
    wb = Workbook()
    for i, r in enumerate(gbd_regions):
        r = r.replace(',', '').replace('/', '_').replace(' ', '_').lower()
        region_keys = []
        for k in keys:
            if k.split(KEY_DELIM_CHAR)[1] == r:
                region_keys.append(k)
        table_region_sheet(dm, region_keys, wb, ("%s.%s" % (i + 1, r)), user, group_size)

    write_data(dm.data, wb)
    write_priors(dm, wb)
        
    f = StringIO()
    wb.save(f)
    f.seek(0)
    return f.read()

def table_region_sheet(dm, keys, wb, name, user, group_size):
    """Make a work sheet for a region

    Parameters
    ----------
    dm_json : str or DiseaseJson object
      the json string or a thin python wrapper around this data that
      is to be plotted
    keys : list the keys to include
    wb : work book
    name : name of the work sheet
    user : username
    group_size : positive integer smaller than 102
    """
    ws = wb.add_sheet(name)
    keys_male_1 = []
    keys_male_2 = []
    keys_female_1 = []
    keys_female_2 = []
    y = keys[0].split(KEY_DELIM_CHAR)[2]
    for k in keys:
        if k.split(KEY_DELIM_CHAR)[3] == 'male':
            if k.split(KEY_DELIM_CHAR)[2] == y:
                keys_male_1.append(k)
            else:
                keys_male_2.append(k)
        else:
            if k.split(KEY_DELIM_CHAR)[2] == y:
                keys_female_1.append(k)
            else:
                keys_female_2.append(k)
    date = strftime("%Y/%m/%d")
    time = strftime("%H:%M:%S")
    ws.write(0, 0, "Dismod III output, date: %s, time: %s, user: %s" % (date, time, user))
    x = 22
    if group_size != 0:
        x = 10 + dismod3.MAX_AGE / group_size

    # tables
    table_disease_model(dm, keys_male_1, ws, 0, 0, group_size)
    table_disease_model(dm, keys_male_2, ws, 0, 39, group_size)
    table_disease_model(dm, keys_female_1, ws, x + 0, 0, group_size)
    table_disease_model(dm, keys_female_2, ws, x + 0, 39, group_size)

    # covariates
    d = 2 * x
    covariates_dict = dm.get_covariates()
    ws.write(2 + d, 0, "Study Level Covariates")
    ws.write(3 + d, 0, "Name")
    ws.write(3 + d, 1, "Rate")
    ws.write(3 + d, 2, "Error")
    ws.write(3 + d, 3, "Value")
    cov_study = covariates_dict['Study_level']
    for i, key in enumerate(cov_study.keys()):
        ws.write(4 + d + i, 0, key)
        check = 'yes'
        if cov_study[key]['rate']['value'] == 0:
            check = 'no'
        ws.write(4 + d + i, 1, check)
        check = 'yes'
        if cov_study[key]['error']['value'] == 0:
            check = 'no'
        ws.write(4 + d + i, 2, check)
        ws.write(4 + d + i, 3, cov_study[key]['value']['value'])
    ws.write(2 + d, 5, "Country Level Covariates")
    ws.write(3 + d, 5, "Name")
    ws.write(3 + d, 6, "Rate")
    ws.write(3 + d, 7, "Error")
    ws.write(3 + d, 8, "Value")
    cov_country = covariates_dict['Country_level']
    for k, key in enumerate(cov_country.keys()):
        ws.write(4 + d + k, 5, key)
        check = 'yes'
        if cov_country[key]['rate']['value'] == 0:
            check = 'no'
        ws.write(4 + d + k, 6, check)
        check = 'yes'
        if cov_country[key]['error']['value'] == 0:
            check = 'no'
        ws.write(4 + d + k, 7, check)
        ws.write(4 + d + k, 8, cov_country[key]['value']['value'])

def table_disease_model(dm, keys, ws, x, y, group_size):
    """Make a table representation of the disease model data and
    estimates provided

    Parameters
    ----------
    dm_json : str or DiseaseJson object
      the json string or a thin python wrapper around this data that
      is to be plotted
    keys : list
      the keys to include
    ws : work sheet
    x : horizontal shift
    y : vertical shift
    group_size : positive integer smaller than 102
    """
    MAX_AGE = dismod3.MAX_AGE
    group_sizes = [1, 4, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10, 16]
    if group_size > 1:
        group_sizes = []
        for i in range(MAX_AGE / group_size):
            group_sizes.append(group_size)
        if MAX_AGE % group_size > 0:
            group_sizes.append(MAX_AGE % group_size)

    data_hash = GBDDataHash(dm.data)
    c = dismod3.utils.KEY_DELIM_CHAR
    type, region, year, sex = keys[0].split(c)

    # add a key: with-condition-death = with-condition-mortality * prevalence * population
    keys.append('with-condition-death' + c + region + c + year + c + sex)
    
    ws.write(x + 2, y, "Condition: %s" % (dm.params['condition']))
    ws.write(x + 3, y, "Region: %s" % (region))
    ws.write(x + 4, y + 1, "%s %s" % (sex.capitalize(), year))
    x += 5
    for i in range(1, 5):
        ws.write(x, y + i, "Data")
    for i in range(5, 9):
        ws.write(x, y + i, "Prior")
    for i in range(9, 37):
        ws.write(x, y + i, "Posterior")
    x += 1
    ws.write(x, y, "Age")
    ws.write(x, y + 1, "Prevalence")
    ws.write(x, y + 2, "Incidence")
    ws.write(x, y + 3, "Remission")
    ws.write(x, y + 4, "Excess Mortality")
    ws.write(x, y + 5, "Prevalence")
    ws.write(x, y + 6, "Incidence")
    ws.write(x, y + 7, "Remission")
    ws.write(x, y + 8, "Excess Mortality")
    ws.write(x, y + 9, "Prevalence")
    ws.write(x, y + 10, "Prevalence")
    ws.write(x, y + 11, "Prevalence")
    ws.write(x, y + 12, "Incidence")
    ws.write(x, y + 13, "Incidence")
    ws.write(x, y + 14, "Incidence")
    ws.write(x, y + 15, "Remission")
    ws.write(x, y + 16, "Remission")
    ws.write(x, y + 17, "Remission")
    ws.write(x, y + 18, "Excess Mortality")
    ws.write(x, y + 19, "Excess Mortality")
    ws.write(x, y + 20, "Excess Mortality")
    ws.write(x, y + 21, "Duration")
    ws.write(x, y + 22, "Duration")
    ws.write(x, y + 23, "Duration")
    ws.write(x, y + 24, "With-condition Mortality")
    ws.write(x, y + 25, "With-condition Mortality")
    ws.write(x, y + 26, "With-condition Mortality")
    ws.write(x, y + 27, "RR Mortality")
    ws.write(x, y + 28, "RR Mortality")
    ws.write(x, y + 29, "RR Mortality")
    ws.write(x, y + 30, "Age of onset")
    ws.write(x, y + 31, "Incidence_x_duration")
    ws.write(x, y + 32, "Incidence_x_duration")
    ws.write(x, y + 33, "Incidence_x_duration")
    ws.write(x, y + 34, "With-condition Death")
    ws.write(x, y + 35, "With-condition Death")
    ws.write(x, y + 36, "With-condition Death")
    x += 1
    ws.write(x, y, "(years)")
    for i in range(1, 10):
        ws.write(x, y + i, "(rate)")
    ws.write(x, y + 10, "lower ui")
    ws.write(x, y + 11, "upper ui")
    ws.write(x, y + 12, "(rate)")
    ws.write(x, y + 13, "lower ui")
    ws.write(x, y + 14, "upper ui")
    ws.write(x, y + 15, "(rate)")
    ws.write(x, y + 16, "lower ui")
    ws.write(x, y + 17, "upper ui")
    ws.write(x, y + 18, "(rate)")
    ws.write(x, y + 19, "lower ui")
    ws.write(x, y + 20, "upper ui")
    ws.write(x, y + 21, "(years)")
    ws.write(x, y + 22, "lower ui")
    ws.write(x, y + 23, "upper ui")
    ws.write(x, y + 24, "(rate)")
    ws.write(x, y + 25, "lower ui")
    ws.write(x, y + 26, "upper ui")
    ws.write(x, y + 27, "(rate)")
    ws.write(x, y + 28, "lower ui")
    ws.write(x, y + 29, "upper ui")
    ws.write(x, y + 30, "(years)")
    ws.write(x, y + 31, "(rate)")
    ws.write(x, y + 32, "lower ui")
    ws.write(x, y + 33, "upper ui")
    ws.write(x, y + 34, "(x1000)")
    ws.write(x, y + 35, "lower ui")
    ws.write(x, y + 36, "upper ui")
    x += 1
    y30 = y + 30
    if group_size == 1:
        for j in range(MAX_AGE):
            ws.write(x + j, y, j)
            ws.write(x + j, y30, j + .5)
    elif group_size == 0:
        start = 0
        end = 0
        for j, s in enumerate(group_sizes):
            start = end
            end = start + s
            if start == 0:
                ws.write(x + j, y, "0")
            elif start == 85:
                ws.write(x + j, y, "85+")
            else:
                ws.write(x + j, y, "%s-%s" % (start, end - 1))
            ws.write(x + j, y30, .5 * (start + end))
    else:
        for j in range(MAX_AGE / group_size + 1):
            start = j * group_size
            end = start + group_size
            if end > MAX_AGE:
                end = MAX_AGE
            ws.write(x + j, y, "%s-%s" % (start, end - 1))
            ws.write(x + j, y30, .5 * (start + end))
    for k in keys:
        type, region, year, sex = k.split(c)
        data_type = clean(type) + ' data'
        data = data_hash.get(data_type, region, year, sex) \
               + data_hash.get(data_type, region, year, 'total')
        column = y
        if type == 'prevalence':
            column = y + 1
        elif type == 'incidence':
            column = y + 2
        elif type == 'remission':
            column = y + 3
        elif type == 'excess-mortality':
            column = y + 4
        else:
            column = -1
        if column != -1:
            data_all = []
            data_weight_all = []
            for j in range(MAX_AGE):
                data_all.append('')
                data_weight_all.append(0)
            for i in range(len(data)):
                start = data[i]['age_start']
                end = data[i]['age_end']
                if end > MAX_AGE:
                    end = MAX_AGE
                for j in range(start, end + 1):
                    p = data[i]['parameter_value'] / float(data[i]['units'])
                    #std = data[i]['standard_error']
                    #age_weight = data[i]['age_weights'][j - start]
                    data_weight = 1
                    #if std != 0:
                        #data_weight = age_weight / std / std
                    #else:
                        #if p != 0:
                            #data_weight = age_weight * 25 / (p**2 * (1 - p)**2)
                    if data_all[j] == '':
                        data_all[j] = p * data_weight
                    else:
                        data_all[j] += p * data_weight
                    data_weight_all[j] += data_weight
            for j in range(MAX_AGE):
                if data_weight_all[j] != 0:
                    data_all[j] = data_all[j] / data_weight_all[j]
            if group_size == 1:
                for j in range(MAX_AGE):
                    ws.write(x + j, column, data_all[j])
            elif group_size == 0:
                start = 0
                end = 0
                for j, gs in enumerate(group_sizes):
                    start = end
                    end = start + gs
                    s = 0
                    n = 0
                    for i in range(start, end):
                        if data_all[i] != '':
                            s += data_all[i]
                            n += 1
                    if n != 0:
                        ws.write(x + j, column, s / n)        
            else:
                for j in range(MAX_AGE / group_size + 1):
                    start = j * group_size
                    end = start + group_size
                    if end > MAX_AGE:
                        end = MAX_AGE
                    s = 0
                    n = 0
                    for i in range(start, end):
                        if data_all[i] != '':
                            s += data_all[i]
                            n += 1
                    if n != 0:
                        ws.write(x + j, column, s / n)
        if type == 'prevalence':
            column = y + 5
        elif type == 'incidence':
            column = y + 6
        elif type == 'remission':
            column = y + 7
        elif type == 'excess-mortality':
            column = y + 8
        else:
            column = -1
        if column != -1:
            if group_size == 1:
                write_table_age_value(dm, k, 'emp_prior_mean', ws, x, column)
            else:
                write_table_group_value(dm, k, 'emp_prior_mean', ws, x, column, group_sizes)
        if type == 'prevalence':
            column = y + 9
        elif type == 'incidence':
            column = y + 12
        elif type == 'remission':
            column = y + 15
        elif type == 'excess-mortality':
            column = y + 18
        elif type == 'duration':
            column = y + 21
        elif type == 'mortality':
            column = y + 24
        elif type == 'relative-risk':
            column = y + 27
        elif type == 'incidence_x_duration':
            column = y + 31
        elif type == 'with-condition-death':
            column = y + 34
        else:
            column = -1
        if column != -1:
            if group_size == 1:
                write_table_age_value(dm, k, 'mean', ws, x, column)
            else:
                write_table_group_value(dm, k, 'mean', ws, x, column, group_sizes)
        if type == 'prevalence':
            column = y + 10
        elif type == 'incidence':
            column = y + 13
        elif type == 'remission':
            column = y + 16
        elif type == 'excess-mortality':
            column = y + 19
        elif type == 'duration':
            column = y + 22
        elif type == 'mortality':
            column = y + 25
        elif type == 'relative-risk':
            column = y + 28
        elif type == 'incidence_x_duration':
            column = y + 32
        elif type == 'with-condition-death':
            column = y + 35
        else:
            column = -1
        if column != -1:
            if group_size == 1:
                write_table_age_value(dm, k, 'lower_ui', ws, x, column)
            else:
                write_table_group_value(dm, k, 'lower_ui', ws, x, column, group_sizes)
        if type == 'prevalence':
            column = y + 11
        elif type == 'incidence':
            column = y + 14
        elif type == 'remission':
            column = y + 17
        elif type == 'excess-mortality':
            column = y + 20
        elif type == 'duration':
            column = y + 23
        elif type == 'mortality':
            column = y + 26
        elif type == 'relative-risk':
            column = y + 29
        elif type == 'incidence_x_duration':
            column = y + 33
        elif type == 'with-condition-death':
            column = y + 36
        else:
            column = -1
        if column != -1:
            if group_size == 1:
                write_table_age_value(dm, k, 'upper_ui', ws, x, column)
            else:
                write_table_group_value(dm, k, 'upper_ui', ws, x, column, group_sizes)

def write_table_age_value(dm, key, item, ws, x, y):
    """Write estimated values into table for all ages

    Parameters
    ----------
    dm : str or DiseaseJson object
      the json string or a thin python wrapper around this data that
      is to be plotted
    key : key
    item : 
    ws : work sheet
    x : horizontal shift
    y : vertical shift
    """
    c = dismod3.utils.KEY_DELIM_CHAR
    type, region, year, sex = key.split(c)
   
    if type == 'with-condition-death':
        key_prevalence = 'prevalence' + c + region + c + year + c + sex
        key_mortality = 'mortality' + c + region + c + year + c + sex
        if len(dm.get_mcmc(item, key_prevalence)) == dismod3.MAX_AGE and \
           len(dm.get_mcmc(item, key_mortality)) == dismod3.MAX_AGE and \
           len(population_by_region_year_sex(region, year, sex)) == dismod3.MAX_AGE:
            for j in range(dismod3.MAX_AGE):
                ws.write(x + j, y, dm.get_mcmc(item, key_prevalence)[j] * \
                         dm.get_mcmc(item, key_mortality)[j] * \
                         population_by_region_year_sex(region, year, sex)[j])
    else:
        if len(dm.get_mcmc(item, key)) == dismod3.MAX_AGE:
            for j in range(dismod3.MAX_AGE):
                ws.write(x + j, y, dm.get_mcmc(item, key)[j])

def write_table_group_value(dm, key, item, ws, x, y, group_sizes):
    """Write estimated values into table for all age_groups

    Parameters
    ----------
    dm : str or DiseaseJson object
      the json string or a thin python wrapper around this data that
      is to be plotted
    key : key
    item : 
    ws : work sheet
    x : horizontal shift
    y : vertical shift
    group_sizes : list of group sizes in order
    """
    c = dismod3.utils.KEY_DELIM_CHAR
    type, region, year, sex = key.split(c)
    start = 0
    end = 0
   
    if type == 'with-condition-death':
        key_prevalence = 'prevalence' + c + region + c + year + c + sex
        key_mortality = 'mortality' + c + region + c + year + c + sex
        if len(dm.get_mcmc(item, key_prevalence)) == dismod3.MAX_AGE and \
           len(dm.get_mcmc(item, key_mortality)) == dismod3.MAX_AGE and \
           len(population_by_region_year_sex(region, year, sex)) == dismod3.MAX_AGE:
            for j, gs in enumerate(group_sizes):
                start = end
                end = start + gs
                raw_rate = dm.get_mcmc(item, key_prevalence)[start:end] * \
                           dm.get_mcmc(item, key_mortality)[start:end] * \
                           population_by_region_year_sex(region, year, sex)[start:end]
                age_indices = []
                for i in range(len(raw_rate)):
                    age_indices.append(i)
                age_weights = dm.get_population(region)[start:end]
                ws.write(x + j, y, rate_for_range(raw_rate, age_indices, age_weights) / gs)
    else:
        if(len(dm.get_mcmc(item, key)) == dismod3.MAX_AGE):
            for j, gs in enumerate(group_sizes):
                start = end
                end = start + gs
                raw_rate = dm.get_mcmc(item, key)[start:end]
                age_indices = []
                for i in range(len(raw_rate)):
                    age_indices.append(i)
                age_weights = dm.get_population(region)[start:end]
                ws.write(x + j, y, rate_for_range(raw_rate, age_indices, age_weights) / gs)
def write_data(data_list, wb):
    """ Write data as a table that can be loaded into dismod"""

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

    additional_keys = sorted(all_keys - set([clean(k) for k in required_keys] + redundant_keys))

    keys = required_keys + additional_keys
    
    for c, k in enumerate(keys):
        if k == 'GBD Region':
            k = 'Region'
        ws.write(0, c, k)
    for r, d in enumerate(sorted(data_list, key=lambda d: d.get('_row'))):
        for c, k in enumerate(keys):
            val = d.get(clean(k), '')
            if val == 'mortality data':
                val = 'with condition mortality data'
            ws.write(r+1, c, val)
            
def write_priors(dm, wb):
    """ Write json for the priors in the workbook, to make results reproducible"""

    ws = wb.add_sheet('priors')
    for r, type in enumerate(['prevalence', 'incidence', 'remission', 'excess-mortality', 'duration', 'relative-risk']):
        prior_str = dm.get_global_priors(type)
        global_priors = dm.get_global_priors(type)
        ws.write(r, 0, type)
        ws.write(r, 1, prior_str)

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
    population = []
    for i in range(dismod3.MAX_AGE):
        population.append(0)
    for iso in countries_for[region]:
        for i in range(dismod3.MAX_AGE):
            if sex == 'all' or sex == 'total':
                population[i] += population_by_age[iso, year, 'male'][i]
                population[i] += population_by_age[iso, year, 'female'][i]
            else:
                population[i] += population_by_age[iso, year, sex][i]
    return population

