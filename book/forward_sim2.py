""" Generate consistent prevalence rates as function of age from
incidence, remission, and excess-mortality"""


### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)


### @export 'initialize-model'

try:
    dm = dismod3.disease_json.DiseaseJson(open('empty_model.json').read())
except IOError:
    dm = dismod3.disease_json.DiseaseJson(open('../empty_model.json').read())

ages = pl.array(dismod3.settings.gbd_ages)
dm.set_param_age_mesh(ages)
dm.vars = dismod3.generic_disease_model.setup(dm)

def set_rate(rate_type, value):
    key = '%s+north_america_high_income+2005+male' % rate_type
    submodel = dm.vars[key]
    submodel['age_coeffs_mesh'].value = \
        pl.log(pl.maximum(dismod3.settings.NEARLY_ZERO, value))

def set_birth_prev(value):
    key = 'bins+north_america_high_income+2005+male'
    submodel = dm.vars[key]
    submodel['initial']['logit_C_0'].value = mc.logit(pl.maximum(dismod3.settings.NEARLY_ZERO, value))


### @export 'congenital'
set_rate('remission', pl.zeros_like(ages))
set_rate('incidence', pl.zeros_like(ages))
set_birth_prev(.2)
set_rate('excess-mortality', 2*(ages/100.)**2)

book_graphics.plot_age_patterns(dm, yticks=[0,.2,.4])
pl.savefig('forward-sim-congenital.pdf')


### @export 'mental'
set_rate('remission', pl.ones_like(ages)*.25)
set_rate('incidence', pl.where((ages > 14) & (ages < 65), (65-ages)*.001, 0.))
set_birth_prev(0)
set_rate('excess-mortality', 1e-4*pl.ones_like(ages))

book_graphics.plot_age_patterns(dm, yticks=[0,.2,.4])
pl.savefig('forward-sim-mental.pdf')


### @export 'old_age'
set_rate('remission', pl.ones_like(ages)*0.)
set_rate('incidence', pl.where(ages > 50, pl.exp((ages-50.)/25.)*.01, 0.))
set_birth_prev(0)
set_rate('excess-mortality', pl.exp(ages/25.)*.01)

book_graphics.plot_age_patterns(dm, yticks=[0,.2,.4])
pl.savefig('forward-sim-old_age.pdf')


### @export 'reproductive'
set_rate('remission', pl.where(ages<45, .1, .35))
set_rate('incidence', pl.where((ages>14) & (ages<65), .025, 0.))
set_birth_prev(0)
set_rate('excess-mortality', pl.where(ages>14, .1, 1e-4))

book_graphics.plot_age_patterns(dm, yticks=[0,.2,.4])
pl.savefig('forward-sim-reproductive.pdf')


### @export 'incidence_pulse'
set_rate('remission', 0*ages)
set_rate('incidence', pl.where(ages==15, .02, 0.))
set_birth_prev(0)
set_rate('excess-mortality', 0*ages)

book_graphics.plot_age_patterns(dm, yticks=[0,.2,.4])
pl.savefig('forward-sim-incidence_pulse.pdf')
