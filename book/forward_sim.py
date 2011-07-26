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


### @export 'initial-rates'

book_graphics.plot_age_patterns(dm)
pl.savefig('initial.png')


### @export 'more-remission'

key = 'remission+north_america_high_income+2005+male'
submodel = dm.vars[key]
submodel['age_coeffs_mesh'].value = \
    pl.log(.15 * pl.ones_like(ages))

book_graphics.plot_age_patterns(dm)
pl.savefig('more-remission.png')


### @export 'increasing-incidence'

key = 'incidence+north_america_high_income+2005+male'
submodel = dm.vars[key]
submodel['age_coeffs_mesh'].value = \
    pl.log(pl.maximum(dismod3.settings.NEARLY_ZERO, .05*ages/100.))

book_graphics.plot_age_patterns(dm)
pl.savefig('increasing-incidence.png')


### @export 'birth_prevalence'

key = 'bins+north_america_high_income+2005+male'
submodel = dm.vars[key]
submodel['initial']['logit_C_0'].value = mc.logit(.2)

book_graphics.plot_age_patterns(dm)
pl.savefig('birth-prevalence.png')


### @export 'smr_2'

key = 'all-cause_mortality+north_america_high_income+2005+male'
m_all = dm.vars[key]

key = 'excess-mortality+north_america_high_income+2005+male'
submodel = dm.vars[key]
submodel['age_coeffs_mesh'].value = \
    pl.log(2*m_all[ages])

book_graphics.plot_age_patterns(dm)
pl.savefig('higher-smr.png')
