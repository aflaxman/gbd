from __future__ import division
import sys
sys.path += ['.', '..', '/homes/peterhm/gbd/', '/homes/peterhm/gbd/book']
import pylab as pl
import pymc as mc
import pandas
import random

import dismod3
reload(dismod3)

def load_new_model(num, area, data_type):
    ''' opens and loads a dismod model number and returns data from Western Europe
    Parameters
    ----------
    num : int
      dismod model number
    area : str
      geographic area that corresponds to dismod3/GBD 2010 hierarchy
    data_type : str
      one of the epidemiologic parameters allowed
      'p', 'i', 'r', 'f', 'pf', 'csmr', 'rr', 'smr', 'X'
    Results
    -------
    model : data.ModelData
      dismod model
    '''
    model = dismod3.data.load('/home/j/Project/dismod/output/dm-%s'%num) 
    model.keep(areas=[area])
    model.input_data = model.get_data(data_type)
    return model

def test_train(model, data_type, replicate):
    ''' splits data into testing and training data sets 
    testing sets have effective sample size = 0 and standard error = inf
    returns the model, the test set indices and the test set indices of the test set
    Parameters
    ----------
    model : data.ModelData
      dismod model
    data_type : str
      one of the epidemiologic parameters allowed
      'p', 'i', 'r', 'f', 'pf', 'csmr', 'rr', 'smr', 'X'
    replicate : int
      integer to be added to the seed
    Results
    -------
    model : data.ModelData
      model with testing set that has effective sample size = 0 and standard error = inf
    test_ix : list
      list of indices that correspond to the test data
    '''
    # save seed
    random.seed(1234567 + replicate)
    # choose random selection of data
    test_ix = random.sample(model.get_data(data_type).index, len(model.get_data(data_type).index)/4)
    # change effective sample size and standard error of random selection
    model.input_data.ix[test_ix, 'effective_sample_size'] = 0
    model.input_data.ix[test_ix, 'standard_error'] = pl.inf
    return model, test_ix

def create_new_vars(model, rate_model, data_type, reference_area, reference_sex, reference_year, iter, thin, burn):
    ''' creates model.vars according to specifications
    Parameters
    ----------
    model : data.ModelData
      dismod model
    rate_model : str
      a rate model
      'neg_binom', 'binom', 'normal', 'log_norm', 'poisson', 'beta'
    data_type : str
      one of the epidemiologic parameters allowed
      'p', 'i', 'r', 'f', 'pf', 'csmr', 'rr', 'smr', 'X'
    reference_area : str
      reference area for model
    reference_sex : str
      reference sex for model
    reference_year : int
      reference year for model
    Results
    -------
    model : data.ModelData
      model with vars built
    '''
    model_vars = dismod3.data.ModelVars()
    model_vars[data_type] = dismod3.data_model.data_model(data_type, model, data_type, 
                                                          reference_area, reference_sex, reference_year, 
                                                          None, None, None, rate_type=rate_model)
    model.vars += model_vars
    model.vars += dismod3.ism.age_specific_rate(model, data_type)
    dismod3.fit.fit_asr(model, data_type, iter=iter, thin=thin, burn=burn)
    return model    

def bias(pred, obs):
    bias = pl.mean(obs - pred)
    return bias

def rmse(pred, obs):
    n = len(obs)
    rmse = pl.sqrt(sum((obs - pred)**2)/n)
    return rmse
    
def mae(pred, obs):
    mae = pl.median(abs(pred - obs))
    return mae
    
def mare(pred, obs):
    mare = pl.median(((pred - obs)/obs)*100)
    return mare
    
def pc(pred_ui, obs):
    wi_ui = list((obs >= pred_ui[:,0]) & (obs <= pred_ui[:,1])).count(1)
    pc = (100*wi_ui)/len(obs)
    return pc