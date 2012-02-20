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

colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f0', '#ffff33']


import age_group_models as agm
reload(agm)

m = {}

f = scipy.interpolate.interp1d([0, 20, 40, 60, 100], [0, .1, .7, .2, 0])
model = agm.simulate_age_group_data(N=25, delta_true=5., pi_true=f)
agm.fit_age_standardizing_model(model)


pl.figure(figsize=(17., 11), dpi=72)

pl.xlabel('Age (Years)', size=32)
pl.ylabel('Rate (Per PY)', size=32)
pl.xticks(size=28)
pl.yticks(size=28)
pl.grid()

pl.plot(model.ages, model.pi_age_true, '-', label='Truth', linewidth=10, color=colors[0], zorder=-2)
pl.axis([-5, 105, 0., 1.])
pl.savefig('/media/windows/t/age_group_truth.png')

graphics.plot_data_bars(model.input_data, 'talk')
pl.axis([-5, 105, 0., 1.])
pl.savefig('/media/windows/t/age_group_data.png')

pl.plot(model.ages, model.vars['mu_age'].stats()['mean'], 'k-', label='Posterior mean', linewidth=10, zorder=10)
pl.axis([-5, 105, 0., 1.])
pl.savefig('/media/windows/t/age_group_standardize_0.png')


pl.plot(model.ages, model.vars['mu_age'].trace().T, '-', color='grey', alpha=.1, linewidth=10, zorder=-1)

pl.plot(model.ages, model.vars['mu_age'].stats()['95% HPD interval'][:,0], 'k-', label='95% HPD interval', linewidth=2, zorder=10)
pl.plot(model.ages, model.vars['mu_age'].stats()['95% HPD interval'][:,1], 'k-', linewidth=2, zorder=10)


pl.legend(fancybox=True, shadow=True, loc='upper right', prop={'size':32})

pl.axis([-5, 105, 0., 1.])
pl.savefig('/media/windows/t/age_group_standardize.png')
pl.show()

