import pylab as pl
import pandas
import sys
sys.path += ['.', '..', '/homes/peterhm/gbd/', '/homes/peterhm/gbd/book'] 

replicates = int(sys.argv[1])

# opens list of all models and rate types used
model_list = pandas.read_csv('/homes/peterhm/gbd/book/model_info.csv')
model_list = pl.array(model_list['model_list'], dtype='i')
rate_types = pandas.read_csv('/homes/peterhm/gbd/book/model_types.csv')
rate_types = pl.array(rate_types['rate_types'])
stats = pandas.read_csv('/homes/peterhm/gbd/book/model_stats.csv')
stats = list(stats['stats'])

# mean results from each model
results = pandas.DataFrame(pl.zeros((len(model_list),(len(stats)-1)*len(rate_types)+1)), index = model_list)
# list of failed models
failures = []

for num,m in enumerate(model_list):
    output = []
    for r,rate in enumerate(rate_types):
        # create column labels for each rate type
        stat_col = []
        for s in stats:
            if s == 'seed': stat_col.append(s)
            else: stat_col.append(s + rate)
        # combine all replicates
        for i in range(replicates):
            try:
                tmp = pandas.read_csv('/clustertmp/dismod/model_comparison_' + str(m) + rate + str(i) + '.csv')
            except IOError:
                tmp = pandas.DataFrame(pl.ones((1,len(stats)))*pl.nan, columns=stat_col)
                fail = pandas.read_csv('/clustertmp/dismod/model_failure_' + str(m) + rate + str(i) + '.csv')
                failures.append(tuple(fail.ix[0,:]))
            # list of dataframes where each dataframe contains the replicates for a rate type
            if i == 0:
                output.append(tmp)
            # and add to that dataframe
            else:
                output[r] = output[r].append(tmp, ignore_index=True)
        # combine all rate_type dataframes of model to be saved
        if r == 0:
            output_save = output[0]
        else:
            output_save = output_save.join(output[r].drop(['seed'],axis=1))
    # save all info for one model
    output_save.to_csv('/homes/peterhm/gbd/book/model_' + str(m) + '.csv')
    # report summary for model
    if num == 0: results.columns = output_save.columns
    results.ix[m,:] = output_save.mean()

# save all results 
results.to_csv('/homes/peterhm/gbd/book/models.csv')
failures = pandas.DataFrame(failures, columns=['model', 'rate_type', 'replicate'])
failures.to_csv('/homes/peterhm/gbd/book/model_failures.csv')