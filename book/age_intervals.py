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

results = {}

### @export 'data'
#dm = dismod3.load_disease_model(15596)  # epilepsy
#dm = dismod3.load_disease_model(16240)  # af

id = 15596
dir = dismod3.settings.JOB_WORKING_DIR % id
fname = '%s/json/dm-%s.json' % (dir, id)
model = data.ModelData.from_gbd_json(fname)

model.input_data['age_end'] += 1  # change year-end to preferred format

### @export 'plot-prevalence-data'
df = model.input_data
df = df[df['data_type'] == 'p'] # select prevalence data
df = df[df['area'] == 'USA']
df = df[df['year_start'] <= 2000]

pl.figure(**book_graphics.half_page_params)
graphics.plot_data_bars(df)
pl.xlabel('Age (Years)')
pl.ylabel('Prevalence (Per 1)')
pl.axis([-2, 102, -.001, .017])
pl.subplots_adjust(left=.1, right=.99, bottom=.15, top=.95)

pl.savefig('epilepsy_ages_intervals.pdf')

### @export 'save-results'

pl.show()
