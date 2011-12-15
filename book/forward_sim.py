""" Generate consistent prevalence rates as function of age from
incidence, remission, and excess-mortality"""


### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc

import consistent_model
import data_simulation

import book_graphics
reload(book_graphics)


### @export 'initialize-model'
def quadratic(a):
    return 1e-6 * (a * (100. - a) + 100.)

def constant(a):
    return .2 * pl.ones_like(a)

val=dict(i=quadratic, f=constant, r=constant)

types = pl.array(['i', 'r', 'f', 'p'])

## generate simulated data
N=0
model = data_simulation.simple_model(N)

for t in types:
    model.parameters[t]['parameter_age_mesh'] = range(0, 101, 20)

sim = consistent_model.consistent_model(model, 'all', 'total', 'all', {})
for t in 'irf':
    for i, k_i in enumerate(sim[t]['knots']):
        sim[t]['gamma'][i].value = pl.log(val[t](k_i))


### @export 'initial-rates'

book_graphics.plot_age_patterns(sim, rate_types='mortality incidence remission prevalence'.split(), yticks=[0,.1,.2])
pl.savefig('initial.pdf')

### @export 'more-remission'

key = 'remission+north_america_high_income+2005+male'
submodel = dm.vars[key]
submodel['age_coeffs_mesh'].value = \
    pl.log(.15 * pl.ones_like(ages))

book_graphics.plot_age_patterns(dm, rate_types='mortality incidence remission prevalence'.split(), yticks=[0,.1,.2])
pl.savefig('more-remission.pdf')


### @export 'increasing-incidence'

key = 'incidence+north_america_high_income+2005+male'
submodel = dm.vars[key]
submodel['age_coeffs_mesh'].value = \
    pl.log(pl.maximum(dismod3.settings.NEARLY_ZERO, .07*(ages/100.)**.5))

book_graphics.plot_age_patterns(dm, rate_types='mortality incidence remission prevalence'.split(), yticks=[0,.2,.4])
pl.savefig('increasing-incidence.pdf')


### @export 'birth_prevalence'

key = 'bins+north_america_high_income+2005+male'
submodel = dm.vars[key]
submodel['initial']['logit_C_0'].value = mc.logit(.2)

book_graphics.plot_age_patterns(dm, rate_types='mortality incidence remission prevalence'.split(), yticks=[0,.4,.8])
pl.savefig('birth-prevalence.pdf')


### @export 'smr_2'

key = 'all-cause_mortality+north_america_high_income+2005+male'
m_all = dm.vars[key]

key = 'excess-mortality+north_america_high_income+2005+male'
submodel = dm.vars[key]
submodel['age_coeffs_mesh'].value = \
    pl.log(10*m_all[ages])

book_graphics.plot_age_patterns(dm, rate_types='mortality incidence remission prevalence'.split(), yticks=[0,.4,.8])
pl.savefig('higher-smr.pdf')
