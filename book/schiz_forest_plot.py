""" Fit data with several rate models and generate forest plot"""

import pylab as pl
import pymc as mc
import simplejson as json

import dismod3
import book_graphics
reload(book_graphics)

pl.figure(**book_graphics.half_page_params)

vars = json.load(open('schiz_forest.json'))
cy = vars['cy']
r = pl.array(vars['r'])
n = pl.array(vars['n'])

sorted_indices = r.argsort().argsort()

se = 1.96*pl.sqrt(r*(1-r)/n)
ms = pl.maximum(3, pl.sqrt(n) / 50.)

for i in range(len(r)):
    pl.errorbar(r[i], sorted_indices[i], xerr=se[i],
                fmt='gs', mew=1, mec='white',
                ms=ms[i])
    pl.text(-.0011, sorted_indices[i], cy[i], ha='right', va='center')
pl.yticks([])

pl.hlines([-1], -1, 1, linewidth=2, linestyle='--', color='k')

results = vars['results']
for i, k in enumerate('binomial beta-binomial poisson negative-binomial normal log-normal offset-log-normal'.split()):
    pi_med = results[k]['quantiles']['50']
    pi_lb = results[k]['95% HPD interval'][0]
    pi_ub = results[k]['95% HPD interval'][1]
    n = pi_med*(1-pi_med) / ((pi_ub - pi_lb)/(2*1.96))**2
    xerr = [
        [pi_med - pi_lb],
        [pi_ub - pi_med]
        ]

    pl.errorbar(pi_med, -(i+2), xerr=xerr,
                fmt='bo', mew=1, mec='white', ms=5)
    pl.text(-.0011, -(i+2), k, ha='right', va='center')

l,r,b,t=pl.axis()
pl.axis([-.001, .014, b-.5, t+.5])
pl.subplots_adjust(left=.2)

pl.savefig('schiz_forest_plot.png')
pl.show()
