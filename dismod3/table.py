""" All methods for reporting Disease Models in tables

Useful High-level Methods::

    table_by_region_year_sex(dm_json, keys, user, group_size)
    table_by_region(dm_json, keys, user, group_size)
    table(dm_json, keys, user, group_size)

Useful Low-level Methods::

    table_region_sheet(dm, keys, wb, name, user, group_size)
    table_disease_model(dm, keys, ws, x, y, group_size)
    write_table_group_value(dm, key, item, ws, x, y, group_size)
"""

import numpy as np
from pyExcelerator import *
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
    wb.save('output.xls')

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
    wb.save('output.xls')

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
    wb.save('output.xls')

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
        elif type == 'case-fatality':
            column = y + 4
        else:
            column = -1
        if column != -1:
            data_all = []
            for j in range(101):
                data_all.append('')
            for i in range(len(data)):
                start = data[i]['age_start']
                end = data[i]['age_end']
                for j in range(start, end + 1):
                    data_all[j] = data[i]['parameter_value']
            if group_size == 1:
                for j in range(101):
                    ws.write(x + j, column, data_all[j])
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
        elif type == 'case-fatality':
            column = y + 8
        else:
            column = -1
        if column != -1 and len(dm.get_mcmc('emp_prior_mean', k)) == 101:
            if group_size == 1:
                for j in range(101):
                    ws.write(x + j, column, dm.get_mcmc('emp_prior_mean', k)[j])
            else:
                write_table_group_value(dm, k, 'emp_prior_mean', ws, x, column, group_size)
        if type == 'prevalence':
            column = y + 9
        elif type == 'incidence':
            column = y + 12
        elif type == 'remission':
            column = y + 15
        elif type == 'case-fatality':
            column = y + 18
        elif type == 'duration':
            column = y + 21
        elif type == 'relative-risk':
            column = y + 27
        else:
            column = -1
        if column != -1:
            if group_size == 1:
                for j in range(0, 101):
                    ws.write(x + j, column, dm.get_mcmc('mean', k)[j])
            else:
                write_table_group_value(dm, k, 'mean', ws, x, column, group_size)
        if type == 'prevalence':
            column = y + 10
        elif type == 'incidence':
            column = y + 13
        elif type == 'remission':
            column = y + 16
        elif type == 'case-fatality':
            column = y + 19
        elif type == 'duration':
            column = y + 22
        elif type == 'relative-risk':
            column = y + 28
        else:
            column = -1
        if column != -1:
            if group_size == 1:
                for j in range(0, 101):
                    ws.write(x + j, column, dm.get_mcmc('lower_ui', k)[j])
            else:
                write_table_group_value(dm, k, 'lower_ui', ws, x, column, group_size)
        if type == 'prevalence':
            column = y + 11
        elif type == 'incidence':
            column = y + 14
        elif type == 'remission':
            column = y + 17
        elif type == 'case-fatality':
            column = y + 20
        elif type == 'duration':
            column = y + 23
        elif type == 'relative-risk':
            column = y + 29
        else:
            column = -1
        if column != -1:
            if group_size == 1:
                for j in range(0, 101):
                    ws.write(x + j, column, dm.get_mcmc('upper_ui', k)[j])
            else:
                write_table_group_value(dm, k, 'upper_ui', ws, x, column, group_size)

def write_table_group_value(dm, key, item, ws, x, y, group_size):
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
    group_size : positive integer smaller than 102
    """     
    region = key.split(dismod3.utils.KEY_DELIM_CHAR)[1]
    for j in range(100 / group_size + 1):
        start = j * group_size
        end = start + group_size
        if end > 101:
            end = 101
        raw_rate = dm.get_mcmc(item, key)[start:end]
        age_indices = []
        for i in range(len(raw_rate)):
            age_indices.append(i)
        age_weights = dm.get_population(region)[start:end]
        ws.write(x + j, y, rate_for_range(raw_rate, age_indices, age_weights))


