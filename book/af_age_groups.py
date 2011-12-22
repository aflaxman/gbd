""" Explore age groups in AF data and Epilipsy"""


import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}

### @export 'data'
dm = dismod3.load_disease_model(15596)  # epilipsy
#dm = dismod3.load_disease_model(16240)  # af

data = dm.filter_data('prevalence+all+all+all')

hist = pl.zeros((101,101))
for d in data:
    hist[d['age_start'], d['age_end']] += 1

most_freq_cnt = hist.max()

### @export 'scatter-prevalence-age-groups'

pl.figure(**book_graphics.half_page_params)
for a_0 in range(101):
    for a_1 in range(101):
        if hist[a_0, a_1] <= 0.:
            continue

        x_i = .5 * (a_0 + a_1)
        y_i = a_1 - a_0 + 1  # epidemiologist notation age 0-4 is five years

        pl.semilogy([x_i], [y_i], 'o', color='none', ms=pl.sqrt(hist[a_0,a_1])*5+2, mec='k', mew=1)
        #pl.text(x_i, y_i, '%d'%hist[a_0, a_1], ha='center', va='center')

for v in [1, 5, 10]:
    pl.plot([-100], [1], 'o', color='none', ms=pl.sqrt(v)*5+2, mec='k', mew=1, label='%d Observations'%v)

pl.legend(loc='upper left', fancybox=True, shadow=True, numpoints=1)

pl.xlabel('Mean of Age Group (Years)')
pl.ylabel('Width of Age Group (Years)')
pl.axis([-5, 110., .6, 500.])
pl.subplots_adjust(left=.1, right=.99, bottom=.15, top=.95)
pl.savefig('af_age_groups_scatter.png')

### @export 'save-results'
book_graphics.save_json('af_age_groups.json', {'most_freq_cnt': most_freq_cnt})
