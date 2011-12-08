import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

faster_run_flag = False

results = {}

import pandas
df = pandas.read_csv('vzv.csv', index_col=None)
df = df[df['Region']=='North Africa/Middle East']
df = df[df['Year End']>=1997]

### @export 'plot-data'
data_bars = zip(df['Age Start'], df['Age End']+1, df['Parameter Value'])

# make lists of x and y points, faster than ploting each bar
# individually
x = []
y = []
for a_0i, a_1i, p_i in data_bars:
    x += [a_0i, a_1i, pl.nan]
    y += [p_i, p_i, pl.nan]

pl.plot(x, y, 'ks-', mew=1, mec='w', ms=4)
pl.axis([-5, 55, 0, 1])
pl.xlabel('Age (Years)')
pl.ylabel('VZV Seroprevalence (Per 1)')
pl.savefig('vzv_data.pdf')


### @export 'forest-plot-age-5'
df5 = df[(df['Age Start'] <= 5) & (df['Age End'] >= 5)]
r = df5['Parameter Value']
n = df5['effective_sample_size']
l = [a.split(',')[0] + ' et al' for a in df5['authors']]
book_graphics.forest_plot(r, n, data_labels=l, xmax=1.0)

pooled = (r*n).sum() / n.sum()
se = pl.sqrt(pooled*(1-pooled)/n.sum())
# TODO: make diamond width of uncertainty
pl.errorbar(pooled, -.5, xerr=[1.96*se], fmt='kd', mew=1, mec='white', ms=15)

pl.vlines([pooled], -.75, 2.25, linewidth=1, linestyle='dashed', color='k')
pl.axis([0, 1, -.75, 2.25])
pl.savefig('vzv_forest.pdf')

