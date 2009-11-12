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
        x = 10 + 100 / group_size
    table_disease_model(dm, keys_male_1, ws, 0, 0, group_size)
    table_disease_model(dm, keys_male_2, ws, 0, 33, group_size)
    table_disease_model(dm, keys_female_1, ws, x, 0, group_size)
    table_disease_model(dm, keys_female_2, ws, x, 33, group_size)

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
    group_sizes = [1, 4, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10, 16]
    if group_size > 1:
        group_sizes = []
        for i in range(101 / group_size):
            group_sizes.append(group_size)
        group_sizes.append(101 % group_size)

    data_hash = GBDDataHash(dm.data)
    type, region, year, sex = keys[0].split(dismod3.utils.KEY_DELIM_CHAR)
    cnt = len(keys)
    ws.write(x + 2, y, "Condition: %s" % (dm.params['condition']))
    ws.write(x + 3, y, "Region: %s" % (region))
    ws.write(x + 4, y + 1, "%s %s" % (sex.capitalize(), year))
    x += 5
    for i in range(1, 5):
        ws.write(x, y + i, "Data")
    for i in range(5, 9):
        ws.write(x, y + i, "Prior")
    for i in range(9, 31):
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
    ws.write(x, y + 24, "With-condition")
    ws.write(x, y + 25, "With-condition")
    ws.write(x, y + 26, "With-condition")
    ws.write(x, y + 27, "RR mortality")
    ws.write(x, y + 28, "RR mortality")
    ws.write(x, y + 29, "RR mortality")
    ws.write(x, y + 30, "Age of onset")
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
    x += 1
    y30 = y + 30
    if group_size == 1:
        for j in range(101):
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
        for j in range(100 / group_size + 1):
            start = j * group_size
            end = start + group_size
            if end > 101:
                end = 101
            ws.write(x + j, y, "%s-%s" % (start, end - 1))
            ws.write(x + j, y30, .5 * (start + end))
    for k in keys:
        type, region, year, sex = k.split(dismod3.utils.KEY_DELIM_CHAR)
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
            for j in range(101):
                data_all.append('')
                data_weight_all.append(0)
            for i in range(len(data)):
                start = data[i]['age_start']
                end = data[i]['age_end']
                for j in range(start, end + 1):
                    p = data[i]['parameter_value']
                    std = data[i]['standard_error']
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
            for j in range(101):
                if data_weight_all[j] != 0:
                    data_all[j] = data_all[j] / data_weight_all[j]
            if group_size == 1:
                for j in range(101):
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
                for j in range(100 / group_size + 1):
                    start = j * group_size
                    end = start + group_size
                    if end > 101:
                        end = 101
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
        if column != -1 and len(dm.get_mcmc('emp_prior_mean', k)) == 101:
            if group_size == 1:
                for j in range(101):
                    ws.write(x + j, column, dm.get_mcmc('emp_prior_mean', k)[j])
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
        elif type == 'relative-risk':
            column = y + 27
        else:
            column = -1
        if column != -1:
            if group_size == 1 and len(dm.get_mcmc('mean', k)) == 101:
                for j in range(0, 101):
                    ws.write(x + j, column, dm.get_mcmc('mean', k)[j])
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
        elif type == 'relative-risk':
            column = y + 28
        else:
            column = -1
        if column != -1 and len(dm.get_mcmc('lower_ui', k)) == 101:
            if group_size == 1:
                for j in range(0, 101):
                    ws.write(x + j, column, dm.get_mcmc('lower_ui', k)[j])
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
        elif type == 'relative-risk':
            column = y + 29
        else:
            column = -1
        if column != -1 and len(dm.get_mcmc('upper_ui', k)) == 101:
            if group_size == 1:
                for j in range(0, 101):
                    ws.write(x + j, column, dm.get_mcmc('upper_ui', k)[j])
            else:
                write_table_group_value(dm, k, 'upper_ui', ws, x, column, group_sizes)

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
    if(len(dm.get_mcmc(item, key)) == 101):
        region = key.split(dismod3.utils.KEY_DELIM_CHAR)[1]
        start = 0
        end = 0
        for j, gs in enumerate(group_sizes):
            start = end
            end = start + gs
            raw_rate = dm.get_mcmc(item, key)[start:end]
            age_indices = []
            for i in range(len(raw_rate)):
                age_indices.append(i)
            age_weights = dm.get_population(region)[start:end]
            ws.write(x + j, y, rate_for_range(raw_rate, age_indices, age_weights))


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
        ws.write(0, c, k)

    for r, d in enumerate(sorted(data_list, key=lambda d: d['_row'])):
        for c, k in enumerate(keys):
            ws.write(r+1, c, d.get(clean(k), ''))
            
def write_priors(dm, wb):
    """ Write json for the priors in the workbook, to make results reproducible"""

    ws = wb.add_sheet('priors')
    for r, type in enumerate(['prevalence', 'incidence', 'remission', 'excess-mortality']):
        prior_str = dm.get_global_priors(type)
        global_priors = dm.get_global_priors(type)
        ws.write(r, 0, type)
        ws.write(r, 1, prior_str)
