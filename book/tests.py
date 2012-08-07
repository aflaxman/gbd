from __future__ import division
import pylab as pl
import pandas

import model_utilities as mu
reload(mu)

# create data 
pred = pandas.DataFrame(pl.arange(10), columns=['mean'])
pred_ui = pandas.DataFrame(pl.hstack((pred-1, pred+1)), columns=['lower','upper'])
obs = pandas.DataFrame(pl.arange(10)+1, columns=['value'])

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