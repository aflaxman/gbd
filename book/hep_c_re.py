""" Estimate prevalence of HCV in North Africa Middle East region, to see effects of priors on country-level variation

"""
# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import colorbrewer
colors = ['#%x%x%x' % col for col in colorbrewer.Set1[3]]

import dismod3
reload(dismod3)


m = {}
for het in 'Slightly Moderately Very'.split():
    model = dismod3.data.ModelData.load(dismod3.settings.JOB_WORKING_DIR % 29772)
    model.keep(areas=['north_africa_middle_east'], sexes=['male', 'total'])
    model.parameters['p']['heterogeneity'] = het
    model.vars += dismod3.ism.age_specific_rate(model, data_type='p',
                                                reference_area='north_africa_middle_east', reference_sex='male', reference_year=1990)
    dismod3.fit.fit_asr(model.vars, iter=20000, burn=5000, thin=15, tune_interval=100) # this should fit all of the asrs available
    model.vars.plot_acorr()
    model.vars.plot_trace()

    m[het] = model


# use results of fit in a plot
pl.figure()
for i, het in enumerate('Slightly Moderately Very'.split()):
    model = m[het]
    alpha = model.vars['p']['alpha']
    y = []
    x = []
    x_lb = []
    x_ub = []
    for j, alpha_i in enumerate(alpha):
        stats = alpha_i.stats()
        y.append(j+.2*i)
        x.append(stats['mean'])
        x_lb.append(stats['mean'] - stats['95% HPD interval'][0])
        x_ub.append(stats['95% HPD interval'][1] - stats['mean'])
    pl.errorbar(x, y, xerr=[x_lb, x_ub], fmt='s', color=colors[i], label='%s ($\\delta=%.2f$)' % (het, pl.exp(model.vars['p']['eta'].stats()['mean'])))

re_names = model.vars['p']['U'].columns
pl.yticks(range(len(re_names)), re_names, rotation=0)
pl.legend(fancybox=True, shadow=True, title='Heterogeneity', numpoints=1)
pl.vlines([0], -.5, len(re_names)+.5)
pl.show()
