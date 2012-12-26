from __future__ import division
import pylab as pl
import pymc as mc
import pandas

import pdb

rate_types = pandas.read_csv('/homes/peterhm/gbd/book/validity/model_types.csv')
rate_types = pl.array(rate_types['rate_types'])
stats = pandas.read_csv('/homes/peterhm/gbd/book/validity/model_stats.csv')
stats = list(stats['stats'])
model_list = pandas.read_csv('/homes/peterhm/gbd/book/validity/model_list.csv')
model_list = list(model_list['model_list'])

def find_best(model_list, stat, rate_types):
    # create output dataframe (model number by rate type)
    output = pandas.DataFrame(pl.zeros((len(model_list),len(rate_types))), index = model_list, columns = rate_types)
    
    # columns of interest in data
    col = []
    for r in rate_types:
        col.append(stat+'_'+r)
    
    for m in model_list:
        data = pandas.read_csv('/homes/peterhm/gbd/book/validity/model_'+ str(m) + '.csv')
        all_col = set(data.columns)
        # delete columns not in use
        data = data.drop(all_col.difference(col), axis=1)
        data = pl.array(data)
        try: 
            if (stat == 'bias') | (stat == 'rmse') | (stat == 'mae') | (stat == 'mare') | (stat == 'time'):
                # find location of minimum in each row
                loc = pl.nanargmin(abs(data),axis=1)
            elif stat == 'pc':
                # find location closest to 95 in each row
                data = abs(data - 95)
                loc = pl.nanargmin(data,axis=1)
            else: 
                print 'not a valid statistic'
                assert 0
        except ValueError:
            loc = []
        output.ix[m,0] = list(loc).count(0)
        output.ix[m,1] = list(loc).count(1)
        output.ix[m,2] = list(loc).count(2)
        output.ix[m,3] = list(loc).count(3)

    return output

# find summary of statistics    
bias_summ = find_best(model_list, 'bias', rate_types)
rmse_summ = find_best(model_list, 'rmse', rate_types)
mae_summ = find_best(model_list, 'mae', rate_types)
mare_summ = find_best(model_list, 'mare', rate_types)
pc_summ = find_best(model_list, 'pc', rate_types)
time_summ = find_best(model_list, 'time', rate_types)

# save results
bias_summ.to_csv('/homes/peterhm/gbd/book/validity/model_summ_bias.csv')
rmse_summ.to_csv('/homes/peterhm/gbd/book/validity/model_summ_rmse.csv')
mae_summ.to_csv('/homes/peterhm/gbd/book/validity/model_summ_mae.csv')
mare_summ.to_csv('/homes/peterhm/gbd/book/validity/model_summ_mare.csv')
pc_summ.to_csv('/homes/peterhm/gbd/book/validity/model_summ_pc.csv')
time_summ.to_csv('/homes/peterhm/gbd/book/validity/model_summ_time.csv')
    

