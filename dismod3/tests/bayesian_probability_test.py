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

# Testing the Bayesian probability code is complicated because (1) the
# model is still being developed and (2) the output is random

# But you've gotta try...

# Instead of using fixtures, I'm going to try creating everything
# here.  I think that this will be clearer, but if it becomes very
# laborious, I should revisit this decision

def create_test_rates(rate_function_str='(age/100.0)**2', rate_type='prevalence data', age_list=range(100), num_subjects=1000):
    import dismod3.models as models

    params = {}
    params['disease'], flag = models.Disease.objects.get_or_create(name='Test Disease')
    params['region'], flag = models.Region.objects.get_or_create(name='World')
    params['rate_type'] = rate_type
    params['sex'] = 'total'
    params['country'] = 'Canada'
    params['epoch_start'] = 2000
    params['epoch_end'] = 2000
    

    rate_list = []

    # TODO: make this safe and robust
    rate_function = eval('lambda age: %s'%rate_function_str)
    for a in age_list:
        params['age_start'] = a
        params['age_end'] = a

        params['denominator'] = num_subjects
        params['numerator'] = rate_function(a) * params['denominator']

        new_rate = models.Rate(**params)
        new_rate.params['Notes'] = 'Simulated data, created using function %s' % rate_function_str
        new_rate.save()
        rate_list.append(new_rate)
        
    return rate_list

def create_test_asrf(rate_function_str='(age/100.0)**2', rate_type='prevalence data', age_list=range(100), num_subjects=1000):
    import dismod3.models as models

    params = {}
    params['disease'], flag = models.Disease.objects.get_or_create(name='Test Disease')
    params['region'], flag = models.Region.objects.get_or_create(name='World')
    params['rate_type'] = rate_type
    params['sex'] = 'total'
    params['notes'] = 'Simulated rate function, created using %s, with age_list %s' % (rate_function_str, str(age_list))

    rf = models.AgeSpecificRateFunction(**params)
    rf.save()
    rf.rates = create_test_rates(rate_function_str, rate_type, age_list, num_subjects)
    rf.save()

    return rf
