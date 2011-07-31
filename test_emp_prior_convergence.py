from pylab import *

import fit_emp_prior
dm = []
try:
    for i in range(5):
        dm.append(fit_emp_prior.fit_emp_prior(15194, 'prevalence'))
except Exception, e:
    print e

import pymc as mc
for dm_i in dm:
    figure(1)
    plot(exp(dm_i.vars['age_coeffs'].stats()['quantiles'][50]), color=cm.spectral(i/10.))

    #mc.Matplot.plot(dm[0].vars['dispersion'], common_scale=False)
    
show()
