import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}

### @export 'data'
# load model
dm = dismod3.load_disease_model(16397)

# set expert priors and other model parameters
dm.params['global_priors']['level_value']['incidence']['age_before'] = 10
dm.params['global_priors']['level_value']['incidence']['age_after'] = 99
#dm.params['global_priors']['smoothness']['incidence']['age_start'] = 10

dm.params['global_priors']['level_bounds']['remission']['upper'] = .05

dm.params['global_priors']['level_value']['prevalence']['age_before'] = 10

dm.params['global_priors']['smoothness']['relative_risk']['amount'] = 'Very'

dm.params['covariates']['Country_level']['LDI_id_Updated_7July2011']['rate']['value'] = 0
dm.params['covariates']['Study_level']['cv_past_year']['rate']['value'] = 1

# clear any fit and priors
dm.clear_fit()
dm.clear_empirical_prior()
dismod3.neg_binom_model.covariate_hash = {}

# initialize model data
#dismod3.neg_binom_model.fit_emp_prior(dm, 'prevalence')
import fit_world
fit_world.fit_world(dm)

import fit_posterior
fit_posterior.fit_posterior(dm, 'north_america_high_income', 'female', '2005')

### @export 'save'
book_graphics.save_json('bipolar.json', vars())
