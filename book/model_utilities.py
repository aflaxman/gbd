import sys
sys.path += ['.', '..'] #['../gbd', '../gbd/book', '../dm3-computation_only/', '../dm3-computation_only/book']
import pylab as pl
import pymc as mc
import pandas
import random

import dismod3
reload(dismod3)

import book_graphics
reload(book_graphics)

def load_new_model(num):
    ''' opens and loads a dismod model number and returns data from Western Europe
    Parameters
    ----------
    num : intc
      dismod model number
    Results
    -------
    model : data.ModelData
      dismod model
    '''
    model = dismod3.data.load('/home/j/Project/dismod/output/dm-%s'%num) 
    model.keep(areas=['europe_western'])
    return model

def tester_trainer(model, data_type):
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
    Results
    -------
    model : data.ModelData
      model with testing set that has effective sample size = 0 and standard error = inf
    tester_ix : list
      list of indices that correspond to the test data
    tester_data_type_ix : list
      list of indices pertaining only to the test data set 
    '''
    random.seed(12345)
    # choose random selection of data
    tester_ix=random.sample(model.get_data(data_type).index, len(model.get_data(data_type).index)/4)
    # find index of just that data type
    model_data_type_ix = list(model.get_data('p').index)
    tester_data_type_ix = []
    for i in tester_ix:
        tester_data_type_ix.append(model_data_type_ix.index(i))
    # change effective sample size and standard error of random selection
    model.input_data.ix[tester_ix, 'effective_sample_size'] = 0
    model.input_data.ix[tester_ix, 'standard_error'] = pl.inf
    return model, tester_ix, tester_data_type_ix

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
