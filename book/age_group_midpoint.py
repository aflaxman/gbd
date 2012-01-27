""" Where midpoint model works and does not
"""

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import pymc as mc
import scipy.interpolate

import graphics
import book_graphics
reload(book_graphics)


import age_group_models as agm
reload(agm)

m = {}

f = scipy.interpolate.interp1d([0, 20, 40, 60, 100], [.4, .425, .6, .5, .4])
model = agm.simulate_age_group_data(N=15, delta_true=5e2, pi_true=f)
agm.fit_midpoint_covariate_model(model)
m[0] = model

f = scipy.interpolate.interp1d([0, 20, 40, 60, 100], [0, .1, .9, .2, 0])
model = agm.simulate_age_group_data(N=15, delta_true=5e2, pi_true=f)
agm.fit_midpoint_covariate_model(model)
m[1] = model


pl.figure(**book_graphics.full_page_params)
for ii in range(2):
    pl.subplot(2, 1, ii+1)
    model = m[ii]
    #
    graphics.plot_data_bars(model.input_data)
    #
    pl.plot(model.ages, model.vars['mu_age'].trace().T, '-', color='grey', alpha=.1)
    pl.plot(model.ages, model.vars['mu_age'].stats()['mean'], 'w-', linewidth=3)
    pl.plot(model.ages, model.vars['mu_age'].stats()['mean'], 'k-', label='Posterior mean')
    pl.plot(model.ages, model.vars['mu_age'].stats()['95% HPD interval'][:,0], 'k:', label='95% HPD interval')
    pl.plot(model.ages, model.vars['mu_age'].stats()['95% HPD interval'][:,1], 'k:')
    #
    pl.plot(model.ages, model.pi_age_true, 'w-', linewidth=3)
    pl.plot(model.ages, model.pi_age_true, 'k--', label='Truth')
    #
    if ii == 0:
        pl.legend(fancybox=True, shadow=True, loc='upper right')
    pl.xlabel('Age (Years)')
    pl.ylabel('Rate (Per PY)')
    pl.axis([-5, 105, 0., 1.])
    # TODO: yticks overlap
    # TODO: label panels (a) and (b)

pl.subplots_adjust(.1, .1, .98, .98, .275, 0)
pl.savefig('age_group_midpoint.pdf')


pl.show()
