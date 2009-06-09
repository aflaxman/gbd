#!/usr/local/bin/python2.5

import dismod3

import sys

dm = dismod3.get_disease_model(894)

if len(sys.argv) == 2:
    k = dm.params['priors'].keys()[int(sys.argv[1])]

    dm.data = [d for d in dm.data if \
               d['gbd_region'] == k.replace('prevalence data+', '')]

    dm.set_priors('prevalence data', ' smooth 10\n zero 0 15\n zero 99 100\n confidence 1000 .0001\n')

    import dismod3.beta_binomial_model as model
    print 'Processing %s (%d data points)' % (k, len(dm.data))
    model.fit(dm, 'map')
    model.fit(dm, 'mcmc')
    model.fit(dm, 'map')
    model.fit(dm, 'mcmc')
    model.fit(dm, 'map')

else:
    keys = dm.params['priors'].keys()
    for k in keys:
        dm.set_priors(k, ' smooth 10\n zero 0 15\n zero 99 100\n confidence 1000 .0001\n')

    import dismod3.multiregion_model as model
    print 'Processing all regions (%d data points)' % len(dm.data)
    model.fit(dm, 'map')
    model.fit(dm, 'mcmc')
    model.fit(dm, 'map')
    model.fit(dm, 'mcmc')
    model.fit(dm, 'map')
    
print dismod3.post_disease_model(dm)


