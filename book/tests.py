from __future__ import division
import pylab as pl
import pandas

import model_utilities as mu
reload(mu)

# create data 
pred = pandas.DataFrame(pl.arange(10), columns=['mean'])
pred_ui = pandas.DataFrame(pl.hstack((pred-1, pred+1)), columns=['lower','upper'])
obs = pandas.DataFrame(pl.arange(10)+1, columns=['value'])
L = [2, 3, 4, 4, 3, 5, 4, 3, 2, 2, 4, 2, 4, 5, 2, 3]

def test_fa1():
    ix = mu.find_all(L, 5)
    assert ix == [5, 13]
    
def test_fa2():
    ix = mu.find_all(L, 4)
    assert ix == [2, 3, 6, 10, 12]

def test_fa3():
    ix = mu.find_all(L, 6)
    assert ix == []    
    
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