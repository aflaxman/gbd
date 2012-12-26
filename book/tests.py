from __future__ import division
import pylab as pl
import pandas
import sys
sys.path += ['.', '..', '/homes/peterhm/gbd/', '/homes/peterhm/gbd/book']

import model_utilities as mu
reload(mu)

import dismod3
reload(dismod3)

model_num = 40418
test_area = 'europe_western'
data_type = 'p'

# example model0, to test vars and test-train
model = mu.load_new_model(model_num, test_area, data_type)
nan_ix = list(model.input_data['effective_sample_size'][pl.isnan(model.input_data['effective_sample_size'])==1].index) # list of nan in effective sample size
model = mu.create_uncertainty(model, 'binom')
for cv in list(model.input_data.filter(like='x_').columns): # fill missing with 0
    model.input_data[cv] = model.input_data[cv].fillna([0])

# example model1, to test test-train
model1 = mu.load_new_model(model_num, test_area, data_type)
model1 = mu.create_uncertainty(model1, 'normal')

# example model2, to test loading and uncertainty
model2 = mu.load_new_model(model_num, test_area, data_type)
nan_ix2 = list(model.input_data['effective_sample_size'][pl.isnan(model.input_data['effective_sample_size'])==1].index) # list of nan in effective sample size
ten_percent = pl.percentile(model2.input_data['effective_sample_size'], 10.)
model2 = mu.create_uncertainty(model2, 'normal')

# find official areas of western europe
areas = [test_area]
for i in model.hierarchy.edges(test_area):
    areas.append(i[1])

# create data for math functions
pred = pandas.DataFrame(pl.arange(10), columns=['mean'])
pred_ui = pandas.DataFrame(pl.hstack((pred-1, pred+1)), columns=['lower','upper'])
obs = pandas.DataFrame(pl.arange(10)+1, columns=['value'])


def test_load_area():
    # find model unique areas
    model_areas = set(pl.unique(model2.input_data['area']))
    # check that only official areas are listed
    assert model_areas.issubset(areas) == 1

def test_load_datatype():
    data_types = list(pl.unique(model2.input_data['data_type']))
    assert data_types == [data_type]

def test_uncertainty_binom():
    assert pl.all(model2.input_data['effective_sample_size'] >= 0)

def test_uncertainty_binom_value():
    assert pl.all(pl.array(model2.input_data.ix[nan_ix,'effective_sample_size']).round(2) == round(ten_percent,2))

def test_uncertainty_normal():
    assert pl.all(model1.input_data['standard_error'] >= 0)

# now that uncertainty tests passed, can test test_train function and creation of vars
model, test_ix = mu.test_train(model, data_type, 23)
model1, test_ix1 = mu.test_train(model1, data_type, 23)
model = mu.create_new_vars(model, 'binom', data_type, test_area, 'male', 2005)

def test_test_train_se():
    assert pl.all(model1.input_data.ix[test_ix1,'standard_error'] == pl.inf)

def test_test_train_ess():
    assert pl.all(model1.input_data.ix[test_ix1,'effective_sample_size'] == 0)

def test_replicate():
    assert test_ix == test_ix1

def test_vars():
    assert model.vars[data_type] != {}

def test_bias():
    bias = mu.bias(pred, obs)
    assert bias == 1

def test_rmse():    
    rmse = mu.rmse(pred, obs)
    assert rmse == 1

def test_mae():
    mae = mu.mae(pred, obs)
    assert mae == 1

def test_mare():
    mare = mu.mare(pred, obs)
    assert mare == (11./60)*100

def test_pc():
    pc = mu.pc(pred_ui, obs)
    assert pc == 100

def test_pc_fail():
    pred_ui.ix[0,'upper'] = 0 # change UI interval from [-1,1] to [-1,0] 
    # now 1 of 10 observations are outside of the predicted UI bounds
    pc = mu.pc(pred_ui, obs)
    assert pc == 90




if __name__ == '__main__':
    import nose
    nose.run()