from __future__ import division
import matplotlib
matplotlib.use('AGG')
import pylab as pl
import pandas
import time
import sys
sys.path += ['.', '..', '/homes/peterhm/gbd/', '/homes/peterhm/gbd/book'] 

import dismod3
reload(dismod3)

import model_utilities as mu
reload(mu)

model_num = int(sys.argv[1])
rate_type = sys.argv[2]
replicate = int(sys.argv[3])

area = 'europe_western'
data_type = 'p'
data_type_full = 'prevalence'

iter=10000
burn=1000
thin=5

# create output structures
stats = ['seed', 'bias_' + rate_type, 'rmse_' + rate_type, 'mae_' + rate_type, 'mare_' + rate_type, 'pc_' + rate_type, 'time_' + rate_type]
output = pandas.DataFrame(pl.zeros((1, len(stats))), columns=stats)
output['seed'] = replicate
failure = []

# load new model
model = mu.load_new_model(model_num, area, data_type)

# fill any missing covariate data with 0s
for cv in list(model.input_data.filter(like='x_').columns):
    model.input_data[cv] = model.input_data[cv].fillna([0])

# replace invalid uncertainty with 10% of data set
model = mu.create_uncertainty(model, rate_type)

# change values of 0 in lognormal model to 1 observation
if rate_type == 'log_normal':
    # find indices where values are 0
    ix = [i for i, x in enumerate(list(model.input_data['value'])) if x == 0]
    # add 1 observation so no values are zero, also change effective sample size
    model.input_data['effective_sample_size'][ix] = model.input_data['effective_sample_size'][ix] + 1
    model.input_data['value'][ix] = 1.0/model.input_data['effective_sample_size'][ix]
elif rate_type == 'log_offset':
    os.chdir('/homes/peterhm/dismod_cpp-20121204/build')
    os.system('bin/get_data.py %s' %model_num)
    os.system('bin/fit.sh %s %s' %(model_num, data_type_full)) # this need 
    os.system('cd fit/%s ../../src/bin/dismod_spline parameter_in.csv prior_in.csv data_in.csv sample_out.csv info_out.csv' %model_num) # default values in input files
    os.system('burn_in=0.1 ../../bin/summary.py $burn_in parameter_in.csv data_in.csv sample_out.csv statistics_out.csv data_out.csv')
    
# withhold 25% of data
model, test_ix = mu.test_train(model, data_type, replicate)

try:
    # create pymc nodes for model and fit the model
    model = mu.create_new_vars(model, rate_type, data_type, area, 'male', 2005)
    # fit the model, using a hill-climbing alg to find an initial value
    # and then sampling from the posterior with MCMC
    start = time.clock()
    dismod3.fit.fit_asr(model, data_type, iter=iter, thin=thin, burn=burn)
    elapsed = (time.clock() - start)
    
    # extract posterior predicted values for data
    pred = pandas.DataFrame(model.vars[data_type]['p_pred'].stats()['mean'], columns=['mean'], index=model.input_data.index)
    pred_ui = pandas.DataFrame(model.vars[data_type]['p_pred'].stats()['95% HPD interval'], columns=['lower', 'upper'], index=model.input_data.index) 
    obs = pandas.DataFrame(model.vars[data_type]['p_obs'].value, columns=['value'], index=model.input_data.index)

    # subset only test data
    pred_test = pred.ix[test_ix]
    obs_test = obs.ix[test_ix]
    pred_ui_test = pred_ui.ix[test_ix]

    # record statistics for test data
    output['bias_'+rate_type] = mu.bias(pred_test, obs_test)
    output['rmse_'+rate_type] = mu.rmse(pred_test, obs_test)
    output['mae_'+rate_type] = mu.mae(pred_test, obs_test)
    output['mare_'+rate_type] = mu.mare(pred_test, obs_test)
    output['pc_'+rate_type] = mu.pc(pred_ui_test, obs_test)
    output['time_'+rate_type] = elapsed
    
    # save information
    output.to_csv('/clustertmp/dismod/model_comparison_' + str(model_num) + rate_type + str(replicate) + '.csv')
    
    # create and save conversion plots
    dismod3.graphics.plot_acorr(model.vars)
    pl.savefig('/clustertmp/dismod/model_comparison_' + str(model_num) + rate_type + str(replicate) + 'acorr.pdf')
    dismod3.graphics.plot_trace(model.vars)
    pl.savefig('/clustertmp/dismod/model_comparison_' + str(model_num) + rate_type + str(replicate) + 'trace.pdf')    

    # save statistic types (only for 1st replicate)
    if replicate == 0:
        model_stats = pandas.DataFrame(['seed', 'bias_', 'rmse_', 'mae_', 'mare_', 'pc_', 'time_'], columns=['stats'])
        model_stats.to_csv('/homes/peterhm/gbd/book/validity/model_stats.csv')
    
except Exception, e:
    print e
    # want to know which models fail 
    failure.append((model_num, rate_type, replicate))
    failure = pandas.DataFrame(failure, columns=['model', 'rate_type', 'replicate'])
    failure.to_csv('/clustertmp/dismod/model_failure_' + str(model_num) + rate_type + str(replicate) + '.csv')
    