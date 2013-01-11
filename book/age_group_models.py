""" Compare approaches to age group modeling
"""

# add to path, to make importing possible
import sys
sys.path += ['..','../tests']
import os
import pylab as pl
import pymc as mc

import book_graphics
reload(book_graphics)
import validate_age_group
reload(validate_age_group)
from validate_age_group import *

# set font
book_graphics.set_font()

fit_age_standardizing_model.fmt = '^-1w7'
fit_midpoint_model.fmt = 'o-1w5'
fit_midpoint_covariate_model.fmt = 'x-2k5'
fit_disaggregation_model.fmt = '*-1w13'


def plot_fits(m):
    # set font
    book_graphics.set_font()
    
    pl.figure(**book_graphics.half_quarter_page_params)
    for ii in range(2):
        pl.subplot(2, 2, ii+1)
        model = m[ii]
        #
        graphics.plot_data_bars(model.input_data, color='grey')
        #
        pl.plot(model.ages, model.pi_age_true, 'w-', linewidth=3)
        pl.plot(model.ages, model.pi_age_true, 'k--', label='Truth')
        #
        #pl.plot(model.ages, model.vars['mu_age'].trace().T, '-', color='grey', alpha=.1)
        pl.plot(model.ages, model.vars['mu_age'].stats()['mean'], 'w-', linewidth=4)
        pl.plot(model.ages, model.vars['mu_age'].stats()['mean'], 'k-', linewidth=2, label='Posterior mean')
        
        #
        if ii == 0:
            pl.legend(fancybox=True, shadow=True, loc='upper center', bbox_to_anchor=(.5,-.33), prop={'size': 'x-large'})
        pl.xlabel('Age (years)')
        pl.ylabel('Rate (per PY)')
        pl.axis([-5, 105, -.05, 1.])
        pl.yticks([0, .25, .5, .75])
        pl.xticks(fontsize='large')
        book_graphics.subtitle('(%s)'%'ab'[ii])

    pl.subplots_adjust(top=.93, bottom=.53, wspace=.35)


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
    pl.xlabel('Age (years)')
    pl.ylabel('Rate (per PY)')
    pl.axis([-5, 105, 0., 1.5])
    pl.subplots_adjust(.1, .175, .98, .875, .275)
    pl.savefig('book/graphics/age_group_models.pdf')

    pl.show()
