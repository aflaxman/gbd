""" Script to compare the results of fitting with heterogeneity slight
and very in the empirical prior phase"""

import pylab as pl

import colorbrewer
colors = ['#%x%x%x' % col for col in colorbrewer.Set1[6]]

import dismod3

m = {}
for het in ['Slightly', 'Moderately', 'Very']:
    model = dismod3.data.fetch_disease_model_if_necessary(30026, 'models/mi/')
    model.keep(areas=['europe_western'], sexes=['male', 'total'], start_year=2000)

    model.parameters['i']['heterogeneity'] = het

    model.vars += dismod3.ism.age_specific_rate(model, 'i')
    dismod3.fit.fit_asr(model, 'i')
    
    m[het] = model

# display uncertainty in age pattern
pl.clf()

for i, het in enumerate(['Slightly', 'Moderately', 'Very']):
    est = m[het].vars['i']['mu_age'].stats()
    #
    x = pl.arange(40,101,10)
    y = est['mean'][x]
    yerr = [y - est['95% HPD interval'][x,0], est['95% HPD interval'][x,1] - y]
    #
    pl.errorbar(x+i-1, y, yerr=yerr, fmt='s-', color=colors[i], mec='w', label=het)

pl.legend(title='Heterogeneity:', fancybox=True, shadow=True, loc='upper left')
pl.xlabel('Age (Years)')
pl.ylabel('Incidence (Per PY)')
pl.title('Parameter Uncertainty of Empirical Prior\n(incidence data for europe_western only)')

dismod3.graphics.plot_data_bars(model.get_data('i'))
pl.axis([35,105,-.001,.11])
