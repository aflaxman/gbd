# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Generate consistent prevalence rates as function of age from incidence, remission, and excess-mortality
# ====================

# <codecell>

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc
import pandas

import consistent_model
import data_simulation

import book_graphics
reload(book_graphics)

# <codecell>

### @export 'initialize-model'
types = pl.array(['i', 'r', 'f', 'p'])
model = data_simulation.simple_model(0)
model.input_data = pandas.read_csv('/home/j/Project/dismod/gbd/data/ssas_mx.csv', index_col=None)

for t in types:
    model.parameters[t]['parameter_age_mesh'] = [0, 5, 15, 25, 35, 45, 55, 65, 75, 100]

model.vars = consistent_model.consistent_model(model, 'all', 'total', 'all', {})
for i, k_i in enumerate(model.parameters[t]['parameter_age_mesh']):
    model.vars['i']['gamma'][i].value = pl.log(k_i*.0001 + .001)
    model.vars['r']['gamma'][i].value = pl.log(.1)
    model.vars['f']['gamma'][i].value = pl.log(.05)

# <codecell>

### @export 'initial-rates'
reload(book_graphics)
book_graphics.plot_age_patterns(model, types='i r m f p'.split(), xticks=[0,50,100],
                                yticks=dict(i=[0,.01,.02], r=[0,.05,.1], m=[0,.2,.4], f=[0,.05,.1], p=[0,.05,.1]))
pl.subplots_adjust(wspace=.5)
pl.savefig('book/graphics/initial.pdf')

# <codecell>

### @export 'more-remission'
reload(book_graphics)
for i, k_i in enumerate(model.parameters[t]['parameter_age_mesh']):
    model.vars['f']['gamma'][i].value = pl.log(k_i*.005 + .001)
book_graphics.plot_age_patterns(model, types='i r m f p'.split(), xticks=[0,50,100],
                                yticks=dict(i=[0,.01,.02], r=[0,.05,.1], m=[0,.2,.4], f=[0,.3,.6], p=[0,.01,.02]),
                                panel='a')
pl.subplots_adjust(wspace=.5)
pl.savefig('book/graphics/more-excess-mortality.pdf')
# <codecell>

### @export 'birth_prevalence'

p_0 = .015
model.vars['logit_C0'].value = mc.logit(p_0)
p = model.vars['p']['mu_age'].value

print """
 For a condition with prevalence of
  %.1f\\%% at age $0$, these rates yield a prevalence age pattern which is
  highly nonlinear, dipping to a minimum of %.1f\\%% at age %d, and then
  increasing back up to %.1f\\%% at the oldest ages.
""" % (p_0*100, p.min()*100, p.argmin(), p[-1]*100)

book_graphics.plot_age_patterns(model, types='i r m f p'.split(), xticks=[0,50,100],
                                yticks=dict(i=[0,.01,.02], r=[0,.05,.1], m=[0,.2,.4], f=[0,.3,.6], p=[.01,.015,.02]),
                                panel='b')
pl.subplots_adjust(wspace=.5)
pl.savefig('book/graphics/birth-prevalence.pdf')
# <codecell>

pl.show()
