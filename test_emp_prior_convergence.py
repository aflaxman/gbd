from pylab import *

import fit_emp_prior
dm = []
try:
    for i in range(5):
        dm.append(fit_emp_prior.fit_emp_prior(15194, 'prevalence'))
except KeyboardInterrupt, e:
    print e

import pymc as mc
for dm_i in dm:
    figure(1)
    plot(exp(dm_i.vars['age_coeffs'].stats()['quantiles'][50]), color=cm.spectral(i/10.))
    plot(exp(dm_i.vars['age_coeffs'].stats()['95% HPD interval']), color=cm.spectral(i/10.))

    mc.Matplot.plot(dm_i.vars['dispersion'], common_scale=False)
    mc.Matplot.plot(dm_i.vars['age_coeffs_mesh'], common_scale=False)
    mc.Matplot.plot(dm_i.vars['study_coeffs'], common_scale=False)
    mc.Matplot.plot(dm_i.vars['region_coeffs'], common_scale=False)

figure(1)
axis([15,100,0,.05])
show()
