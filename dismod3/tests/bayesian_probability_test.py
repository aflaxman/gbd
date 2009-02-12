"""
I'm putting doctests here until I figure out how to run them in the models
>>> probabilistic_utils.const_func([1,2,3], 17.0)
array([ 17., 17., 17.])

>>> M,C = probabilistic_utils.uninformative_prior_gp()
>>> M([0,10,100])
array([-10., -10., -10.])

>>> probabilistic_utils.gp_interpolate([0,100], [0,100], [0,100])
array([ 0., 100.])
"""
