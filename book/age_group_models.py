""" Compare approaches to age group modeling
"""

# add to path, to make importing possible
import sys
sys.path += ['.', '..', '../tests']

import pylab as pl
import pymc as mc

import book_graphics
from validate_age_group import *

fit_age_standardizing_model.fmt = '^-1w7'
fit_midpoint_model.fmt = 'o-1w5'
fit_midpoint_covariate_model.fmt = 'x-2k5'
fit_disaggregation_model.fmt = '*-1w13'


def plot_fits(m):
    pl.figure(**book_graphics.full_page_params)
    for ii in range(2):
        pl.subplot(2, 1, ii+1)
        model = m[ii]
        #
        graphics.plot_data_bars(model.input_data)
        #
        pl.plot(model.ages, model.pi_age_true, 'w-', linewidth=3)
        pl.plot(model.ages, model.pi_age_true, 'k--', label='Truth')
        #
        #pl.plot(model.ages, model.vars['mu_age'].trace().T, '-', color='grey', alpha=.1)
        pl.plot(model.ages, model.vars['mu_age'].stats()['mean'], 'w-', linewidth=4)
        pl.plot(model.ages, model.vars['mu_age'].stats()['mean'], 'k-', linewidth=2, label='Posterior mean')
        pl.plot(model.ages, model.vars['mu_age'].stats()['95% HPD interval'][:,0], 'k-', label='95% HPD interval')
        pl.plot(model.ages, model.vars['mu_age'].stats()['95% HPD interval'][:,1], 'k-')
        #
        if ii == 0:
            pl.legend(fancybox=True, shadow=True, loc='upper right', prop={'size': 'x-large'})
        pl.xlabel('Age (Years)', fontsize='x-large')
        pl.ylabel('Rate (Per PY)', fontsize='x-large')
        pl.axis([-5, 105, -.05, 1.])
        pl.yticks([0, .25, .5, .75], fontsize='large')
        pl.xticks(fontsize='large')
        pl.text(0, .9, '(%s)'%'ab'[ii], va='top', fontsize='x-large')

    pl.subplots_adjust(.1, .1, .98, .98, .275, 0)


if __name__ == '__main__':
    m = {}
    for fit in [fit_midpoint_covariate_model, fit_age_standardizing_model, fit_midpoint_model, fit_disaggregation_model]:
        mc.np.random.seed(1234567)
        m[fit] = simulate_age_group_data(N=30, delta_true=5)
        fit(m[fit])

        model = m[fit]

    pl.figure(**book_graphics.half_page_params)
    graphics.plot_data_bars(model.input_data)
    pl.plot(model.ages, model.pi_age_true, 'w-', linewidth=3)
    pl.plot(model.ages, model.pi_age_true, 'k--', label='Truth')
    for fit in [fit_age_standardizing_model, fit_midpoint_model, fit_midpoint_covariate_model, fit_disaggregation_model]:
        pl.plot(model.ages[::10], m[fit].vars['mu_age'].stats()['mean'][::10], 'w-', linewidth=3)
        pl.plot(model.ages[::10], m[fit].vars['mu_age'].stats()['mean'][::10], marker=fit.fmt[0], linestyle=fit.fmt[1],
                mew=float(fit.fmt[2]), mec=fit.fmt[3], ms=float(fit.fmt[4:]), color='k', label=fit.__doc__)
    pl.legend(fancybox=True, shadow=True, loc='upper right')
    pl.xlabel('Age (Years)')
    pl.ylabel('Rate (Per PY)')
    pl.axis([-5, 105, 0., 1.5])
    pl.subplots_adjust(.1, .175, .98, .875, .275)
    pl.savefig('age_group_models.pdf')

    pl.show()
