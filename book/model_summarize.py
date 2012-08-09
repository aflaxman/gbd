from __future__ import division
import pylab as pl
import pymc as mc
import pandas

rate_types = pandas.read_csv('/homes/peterhm/gbd/book/model_types.csv')
rate_types = pl.array(rate_types['rate_types'])
stats = pandas.read_csv('/homes/peterhm/gbd/book/model_stats.csv')
stats = list(stats['stats'])

data = pandas.read_csv('/homes/peterhm/gbd/book/models.csv')
model_list = list(data.index)

def find_best(model_list, stat, rate_types):
    # create output dataframe (model number by rate type)
    output = pandas.DataFrame(pl.zeros((len(model_list),len(rate_types))), index = model_list, columns = rate_types)
    
    # columns of interest in data
    col = []
    for r in rate_types:
        col.append(stat+'_'+r)
    
    for m in model_list:
        data = pandas.read_csv('/homes/peterhm/gbd/book/model_'+ str(m) + '.csv')
        all_col = set(data.columns)
        # delete columns not in use
        data = data.drop(all_col.difference(col), axis=1)
        if (stat == 'bias') | (stat == 'rmse') | (stat == 'mae') | (stat == 'mare') | (stat == 'time'):
            # find location of minimum in each row
            loc = list(abs(pl.array(data)).argmin(axis=1))
        elif stat == 'pc':
            # find location of maximum in each row
            loc = list(abs(pl.array(data)).argmax(axis=1))
        else: 
            print 'not a valid statistic'
            assert 0
        output.ix[m,0] = loc.count(0)
        output.ix[m,1] = loc.count(1)
        output.ix[m,2] = loc.count(2)
        output.ix[m,3] = loc.count(3)

    return output
    
bias_summ = find_best(model_list, 'bias', rate_types)
rmse_summ = find_best(model_list, 'rmse', rate_types)
mae_summ = find_best(model_list, 'mae', rate_types)
mare_summ = find_best(model_list, 'mare', rate_types)
pc_summ = find_best(model_list, 'pc', rate_types)
time_summ = find_best(model_list, 'time', rate_types)
    


