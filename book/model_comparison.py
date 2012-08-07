from __future__ import division
import pylab as pl
import pandas
import sys
sys.path += ['.', '..', '/homes/peterhm/gbd/', '/homes/peterhm/gbd/book'] 

import dismod3
reload(dismod3)

import model_utilities as mu
reload(mu)

rate_types = ['neg_binom', 'normal', 'log_normal', 'binom']
stats = ['seed']
for r in rate_types:
    stats.append('bias_' + r)
    stats.append('rmse_' + r)
    stats.append('mae_' + r)
    stats.append('mare_' + r)
    stats.append('pc_' + r)

draws = int(sys.argv[1])
model_num = int(sys.argv[2])

area = 'europe_western'
data_type = 'p'

iter=10000
burn=1000
thin=5

output = pandas.DataFrame(pl.zeros((len(range(draws)), len(stats))), columns=stats)

for i in range(draws):
    # load new model
    model = mu.load_new_model(model_num, area, data_type)
    # withhold 25% of data, save seed
    model, test_ix = mu.test_train(model, data_type, i)
    output.ix[i, 'seed'] = i
    for r in rate_types:
        # create pymc nodes for model and fit the model
        model = mu.create_new_vars_fit(model, r, data_type, area, 'male', 2005, iter, thin, burn)

        # extract posterior predicted values for data
        pred = pandas.DataFrame(model.vars[data_type]['p_pred'].stats()['mean'], columns=['mean'], index=model.input_data.index)
        pred_ui = pandas.DataFrame(model.vars[data_type]['p_pred'].stats()['95% HPD interval'], columns=['lower', 'upper'], index=model.input_data.index) 
        obs = pandas.DataFrame(model.vars[data_type]['p_obs'].value, columns=['value'], index=model.input_data.index)
        
        # subset only test data
        pred_test = pred.ix[test_ix]
        obs_test = obs.ix[test_ix]
        pred_ui_test = pred_ui.ix[test_ix]
        
        # methods of comparison
        output.ix[i, 'bias_'+r] = mu.bias(pred_test, obs_test)
        output.ix[i, 'rmse_'+r] = mu.rmse(pred_test, obs_test)
        output.ix[i, 'mae_'+r] = mu.mae(pred_test, obs_test)
        output.ix[i, 'mare_'+r] = mu.mare(pred_test, obs_test)
        output.ix[i, 'pc_'+r] = mu.pc(pred_test, obs_test)

output.to_csv('model_comparison_' + str(model_num) + '.csv')   
