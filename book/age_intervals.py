""" Explore age groups in AF data and Epilipsy"""


import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import graphics
import data
reload(graphics)

import book_graphics
reload(book_graphics)

# set font
book_graphics.set_font()

results = {}

### @export 'data'
#model = data.load('models/af')
model = data.load('/home/j/Project/dismod/output/dm-39544')

model.input_data['age_end'] += 1  # change year-end to preferred format

### @export 'plot-prevalence-data'
df = model.input_data
df = df[df['data_type'] == 'p'] # select prevalence data
df = df[df['area'] == 'USA']
#df = df[df['year_start'] <= 2000]

pl.figure(**book_graphics.half_page_params)
pl.subplot(1,2,2)
graphics.plot_data_bars(df)
pl.xticks()
pl.yticks()
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (per 1)',rotation=0)
pl.axis([-2, 102, -.01, .22])
pl.subplots_adjust(left=.1, right=.99, bottom=.15, top=.95)

pl.savefig('book/graphics/af_ages_intervals.pdf')

print 'The systematic review of the descriptive epidemiology of atrial fibrilation included $%d$ observations of disease prevalence for the United States' % len(df)

### @export 'save-results'

pl.show()
