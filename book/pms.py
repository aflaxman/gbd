import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}

### @export 'data'
dm = dismod3.load_disease_model(16314)
dm.calc_effective_sample_size(dm.data)
some_data = ([d for d in dm.data
              if d['data_type'] == 'prevalence data'
              and d['sex'] == 'female'
              and d['effective_sample_size'] > 1])

# TODO: replace fake year data with real year data (called Year Start (original) and End)

countries = pl.unique([s['region'] for s in some_data])
min_year = min([s['year_start'] for s in some_data])
max_year = max([s['year_end'] for s in some_data])
cy = ['%s-%d'%(s['region'], s['year_start']) for s in some_data]

n = pl.array([s['effective_sample_size'] for s in some_data])
r = pl.array([dm.value_per_1(s) for s in some_data])

s = pl.sqrt(r * (1-r) / n)

min_rate_per_100 = '%d' % round(min([dm.value_per_1(d) for d in some_data])*100)

### @export 'empirical-priors'
dm.vars = dismod3.generic_disease_model(dm, '%s+world+2005+female')


### @export 'save'
book_graphics.save_json('pms.json', vars())
