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
    model = mu.load_new_model(model_num, area, data_type)
    model, test_ix = mu.test_train(model, data_type, i)
    output.ix[i, 'seed'] = i
    for r in rate_types:
        # TODO: update these comments
        # create pymc nodes for model
        model = mu.create_new_vars(model, r, data_type, 'europe_western', 'male', 2005, iter, thin, burn)

        # fit the model, using a hill-climbing alg to find an initial value
        # and then sampling from the posterior with MCMC
        dismod3.fit.fit_asr(model, data_type, iter=iter, thin=thin, burn=burn)
        
        # extract posterior predicted values for data
        pred = model.vars[data_type]['p_pred'].stats()['mean'] 
        pred_ui = model.vars[data_type]['p_pred'].stats()['95% HPD interval']
        obs = model.vars[data_type]['p_obs'].value
        n = model.vars[data_type]['p_pred'].stats()['n']
        
        # testing
        
        output.ix[i, 'bias_'+r] = pl.mean(obs[test_ix] - pred[test_ix])
        output.ix[i, 'rmse_'+r] = pl.sqrt(sum((obs[test_ix] - pred[test_ix])**2)/n)
        output.ix[i, 'mae_'+r] = pl.median(abs(pred[test_ix]-obs[test_ix]))
        output.ix[i, 'mare_'+r] = pl.median(((pred[test_ix]-obs[test_ix])/obs[tester_data_type_ix])*100)
        output.ix[i, 'pc_'+r] = (100*list((pred[test_ix] >= pred_ui[test_ix,0]) & (pred[test_ix] <= pred_ui[test_ix,1])).count(1))/len(test_ix)

output.to_csv('model_comparison_' + str(model_num) + '.csv')   
