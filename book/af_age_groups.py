""" Explore age groups in AF data"""


import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}

### @export 'data'
dm = dismod3.load_disease_model(15596)

data = dm.filter_data('prevalence+all+all+all')

hist = pl.zeros((101,101))
for d in data:
    hist[d['age_start'], d['age_end']] += 1

most_freq_cnt = hist.max()

### @export 'scatter-prevalence-age-groups'
x = []
y = []
for d in data:
    x_i = .5 * (d['age_start'] + d['age_end'])
    y_i = d['age_end'] - d['age_start'] + 1  # epidemiologist notation age 0-4 is five years
    if y_i < 1:
        import pdb; pdb.set_trace()
    x.append(x_i+pl.exp(pl.rand()*.5))
    y.append(y_i+pl.rand()*.5)

pl.figure(**book_graphics.quarter_page_params)
pl.semilogy(x, y, 'gs',
            alpha=1.,
            mec='white', ms=5)

pl.savefig('af_age_groups_scatter.png')

pl.xlabel('Middle of Age Group (Years)')
pl.ylabel('Spread of Age Group (Years)')
pl.axis([.9, 110., .9, 110.])
pl.subplots_adjust(left=.1, right=.99, bottom=.15, top=.95)
pl.savefig('af_age_groups_scatter.png')

### @export 'save-results'
book_graphics.save_json('af_age_groups.json', vars())
