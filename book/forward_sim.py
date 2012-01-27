""" Generate consistent prevalence rates as function of age from
incidence, remission, and excess-mortality"""


### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc
import pandas

import consistent_model
import data_simulation

import book_graphics
reload(book_graphics)


### @export 'initialize-model'
types = pl.array(['i', 'r', 'f', 'p'])
model = data_simulation.simple_model(0)
model.input_data = pandas.read_csv('ssas_mx.csv')

for t in types:
    model.parameters[t]['parameter_age_mesh'] = [0, 5, 15, 25, 35, 45, 55, 65, 75, 100]

model.vars = consistent_model.consistent_model(model, 'all', 'total', 'all', {})
for i, k_i in enumerate(model.parameters[t]['parameter_age_mesh']):
    model.vars['i']['gamma'][i].value = pl.log(k_i*.0001 + .001)
    model.vars['r']['gamma'][i].value = pl.log(.1)
    model.vars['f']['gamma'][i].value = pl.log(.05)

### @export 'initial-rates'

book_graphics.plot_age_patterns(model, types='i r m f p'.split(),
                                yticks=dict(i=[0,.01,.02], r=[0,.05,.1], m=[0,.2,.4], f=[0,.05,.1], p=[0,.05,.1]))
pl.savefig('initial.pdf')

### @export 'more-remission'

for i, k_i in enumerate(model.parameters[t]['parameter_age_mesh']):
    model.vars['f']['gamma'][i].value = pl.log(k_i*.005 + .001)
book_graphics.plot_age_patterns(model, types='i r m f p'.split(),
                                yticks=dict(i=[0,.01,.02], r=[0,.05,.1], m=[0,.2,.4], f=[0,.3,.6], p=[0,.01,.02]))
pl.savefig('more-excess-mortality.pdf')

### @export 'birth_prevalence'

model.vars['logit_C0'].value = mc.logit(.02)
book_graphics.plot_age_patterns(model, types='i r m f p'.split(),
                                yticks=dict(i=[0,.01,.02], r=[0,.05,.1], m=[0,.2,.4], f=[0,.3,.6], p=[.01,.015,.02]))
pl.savefig('birth-prevalence.pdf')

