""" Fit data with several rate models and generate forest plot"""

import pylab as pl
import pymc as mc
import simplejson as json

import dismod3
import book_graphics
reload(book_graphics)

pl.figure(**book_graphics.half_page_params)

l = .25
w = .6
data_b = .5
h = .9

ms_scale=25

vars = json.load(open('schiz_forest.json'))
r = pl.array(vars['r'])
n = pl.array(vars['n'])

sorted_indices = r.argsort().argsort()

se = 1.96*pl.sqrt(r*(1-r)/n)
ms = pl.sqrt(n) / ms_scale

#ax = pl.axes([l, data_b, w, .9-data_b], frameon=False)
for i in range(len(r)):
    pl.errorbar(r[i], sorted_indices[i], xerr=se[i],
                fmt='gs', mew=1, mec='white',
                ms=ms[i])
    pl.text(-.001, sorted_indices[i], 'row %d '%i, ha='right', va='center')
#pl.xticks([])
pl.yticks([])

results = vars['results']
for i, k in enumerate('binomial beta-binomial poisson negative-binomial normal log-normal offset-log-normal'.split()):
    #pl.axes([l, data_b-(i+1)*(h-data_b)/len(results), w, (h-data_b)/len(results)],
    #        frameon=False,
    #        sharex=ax)
    #pl.hist(results[k]['trace'], bins=50, normed=True, histtype='stepfilled')
    #pl.xticks([])
    #pl.yticks([])
    pi_med = results[k]['quantiles']['50']
    pi_lb = results[k]['95% HPD interval'][0]
    pi_ub = results[k]['95% HPD interval'][1]
    n = pi_med*(1-pi_med) / ((pi_ub - pi_lb)/(2*1.96))**2
    ms = pl.sqrt(n)/ms_scale
    xerr = [
        [pi_med - pi_lb],
        [pi_ub - pi_med]
        ]
    pl.errorbar(pi_med, -(i+1.5), xerr=xerr,
                fmt='bo', mew=1, mec='white', ms=5)
    pl.text(-.001, -(i+1.5), k + ' ', ha='right', va='center')

l,r,b,t=pl.axis()
pl.axis([0., r, b-.5, t+.5])
pl.subplots_adjust(left=.2)

pl.savefig('schiz_forest_plot.png')
pl.show()
