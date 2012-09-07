from __future__ import division
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
# withhold 25% of data
model, test_ix = mu.test_train(model, data_type, replicate)

try:
    fail = 0
    # create pymc nodes for model and fit the model
    model = mu.create_new_vars(model, rate_type, data_type, area, 'male', 2005, iter, thin, burn)
    # fit the model, using a hill-climbing alg to find an initial value
    # and then sampling from the posterior with MCMC
    start = time.clock()
    dismod3.fit.fit_asr(model, data_type, iter=iter, thin=thin, burn=burn)
    elapsed = (time.clock() - start)
except:
    print 'exception'
    try:
        print 'trying again'
        # load new model (if value ==0, want to change value and effective sample size)
        model = mu.load_new_model(model_num, area, data_type)
        
        # find indices where values are 0
        ix = mu.find_all(list(model.input_data['value']), 0)
        # add 1 observation so no values are zero, also change effective sample size
        model.input_data['effective_sample_size'][ix] = model.input_data['effective_sample_size'][ix] + 1
        model.input_data['value'][ix] = 1.0/model.input_data['effective_sample_size'][ix]
        
        # withhold 25% of data
        model, test_ix = mu.test_train(model, data_type, replicate)
        
        # create pymc nodes for model and fit the model
        model = mu.create_new_vars(model, rate_type, data_type, area, 'male', 2005, iter, thin, burn)
        # fit the model, using a hill-climbing alg to find an initial value
        # and then sampling from the posterior with MCMC
        start = time.clock()
        dismod3.fit.fit_asr(model, data_type, iter=iter, thin=thin, burn=burn)
        elapsed = (time.clock() - start)
    except:
        print 'fail'
        fail = 1
        # want to know which models fail 
        failure.append((model_num, rate_type, replicate))
        failure = pandas.DataFrame(failure, columns=['model', 'rate_type', 'replicate'])
        failure.to_csv('/clustertmp/dismod/model_failure_' + str(model_num) + rate_type + str(replicate) + '.csv')

if fail != 1:      
    # extract posterior predicted values for data
    pred = pandas.DataFrame(model.vars[data_type]['p_pred'].stats()['mean'], columns=['mean'], index=model.input_data.index)
    pred_ui = pandas.DataFrame(model.vars[data_type]['p_pred'].stats()['95% HPD interval'], columns=['lower', 'upper'], index=model.input_data.index) 
    obs = pandas.DataFrame(model.vars[data_type]['p_obs'].value, columns=['value'], index=model.input_data.index)

    # subset only test data
    pred_test = pred.ix[test_ix]
    obs_test = obs.ix[test_ix]
    pred_ui_test = pred_ui.ix[test_ix]

    # methods of comparison
    output['bias_'+rate_type] = mu.bias(pred_test, obs_test)
    output['rmse_'+rate_type] = mu.rmse(pred_test, obs_test)
    output['mae_'+rate_type] = mu.mae(pred_test, obs_test)
    output['mare_'+rate_type] = mu.mare(pred_test, obs_test)
    output['pc_'+rate_type] = mu.pc(pred_ui_test, obs_test)
    output['time_'+rate_type] = elapsed
    
    # save information
    output.to_csv('/clustertmp/dismod/model_comparison_' + str(model_num) + rate_type + str(replicate) + '.csv')
    
    # create and save conversion plots
    model.vars.plot_acorr()
    savefig('/clustertmp/dismod/model_comparison_' + str(model_num) + rate_type + str(replicate) + 'acorr.pdf')
    model.vars.plot_trace()
    savefig('/clustertmp/dismod/model_comparison_' + str(model_num) + rate_type + str(replicate) + 'trace.pdf')    

    # save statistic types (only for 1st replicate)
    if replicate == 0:
        model_stats = pandas.DataFrame(['seed', 'bias_', 'rmse_', 'mae_', 'mare_', 'pc_', 'time_'], columns=['stats'])
        model_stats.to_csv('/homes/peterhm/gbd/book/validity/model_stats.csv')