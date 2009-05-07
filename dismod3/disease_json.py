import numpy as np
import pymc as mc
from pymc import gp

import twill.commands as twc
import simplejson as json

from dismod3.settings import *

from bayesian_models.probabilistic_utils import trim, uninformative_prior_gp, NEARLY_ZERO, MAX_AGE
MISSING = -99


class DiseaseJson:
    def __init__(self, json_str):
        dm = json.loads(json_str)
        self.params = dm['params']
        self.data = dm['data']

    def to_json(self):
        return json.dumps({'params': self.params, 'data': self.data})

    def set_key_by_type(self, key, type, value):
        if not self.params.has_key(key):
            self.params[key] = {}
        self.params[key][type] = value
    def get_key_by_type(self, key, type, default=None):
        return self.params.get(key, {}).get(type, default)

    def get_initial_value(self, type):
        return np.array(self.get_key_by_type('normal_approx', type))
    def set_initial_value(self, type, val):
        self.set_key_by_type('normal_approx', type, list(val))

    def get_map(self, type):
        return np.array(self.get_key_by_type('map', type, default=[]))
    def set_map(self, type, val):
        self.set_key_by_type('map', type, list(val))

    def get_mcmc(self, est_type, data_type):
        return np.array(self.get_key_by_type('mcmc_%s' % est_type, data_type, default=[]))
    def set_mcmc(self, est_type, data_type, val):
        self.set_key_by_type('mcmc_%s' % est_type, data_type, list(val))

    def get_units(self, type):
        return self.get_key_by_type('units', type)
    def set_units(self, type, units):
        self.set_key_by_type('units', type, units)

    def get_priors(self, type):
        return self.get_key_by_type('priors', type) or ''
    def set_priors(self, type, priors):
        self.set_key_by_type('priors', type, priors)

    def get_estimate_age_mesh(self):
        return self.params.get('estimate_age_mesh', [0, 100])
    def set_estimate_age_mesh(self, mesh):
        self.params['estimate_age_mesh'] = list(mesh)

    def get_param_age_mesh(self):
        return self.params['param_age_mesh']
    def set_param_age_mesh(self, mesh):
        self.params['param_age_mesh'] = list(mesh)

    def set_model_source(self, source_obj):
        try:
            import inspect
            self.params['model_source'] = inspect.getsource(source_obj)
        except:
            self.params['model_source'] = '(failed to read file)'
    def get_model_source(self):
        return self.params.get('model_source', '')
    
    def filter_data(self, data_type=None, sex=None):
        return [d for d in self.data if ((not data_type) or d['data_type'] == data_type) \
                                    and ((not sex) or d['sex'] == sex)
                ]
    def extract_units(self, d):
        """
        d is a data hash which might include
        the key 'units', which is a decription
        of the units for this datum.
        
        return the float that d['value'] should
        be multiplied to make the units per 1.0
        """
        try:
            unit_str = d.get('units', '1')
            if unit_str.strip()[0:4] == 'per ':
                units = 1. / float(unit_str.split()[-1])
            else:
                units = float(unit_str)
            return units
        except ValueError:
            print 'could not parse unit str: %s' % unit_str
            return 1.


    def mortality(self):
        """
        calculate the all-cause mortality rate for the
        region and sex of disease_model, and return it
        in an array corresponding to age_mesh
        """
        mortality_data = self.filter_data('all-cause mortality data')
        if len(mortality_data) == 0:
            return np.zeros(len(self.get_estimate_age_mesh()))
        else:
            self.fit_normal_approx('all-cause mortality data')
            return self.get_initial_value('all-cause mortality data')

    def value_per_1(self, data):
        scale = self.extract_units(data)
        return data['value'] * scale

    def se_per_1(self, data):
        scale = self.extract_units(data)
        if data['standard_error'] == MISSING:
            return MISSING
        else:
            se = data['standard_error']
            
        se *= scale

        return se
        

    def fit_normal_approx(self, data_type):
        """
        This 'normal approximation' estimate for an age-specific dataset
        is formed by using each datum to produce an estimate of the
        function value at a single age, and then saying that the logit of
        the true rate function is a gaussian process and these
        single age estimates are observations of this gaussian process.
        
        This is less valid and less accurate than using MCMC or MAP, but
        it is much faster.  It is used to generate an initial value for
        the maximum-liklihood estimate.
        """
        from bayesian_models import probabilistic_utils
        from bayesian_models.probabilistic_utils import trim, uninformative_prior_gp, NEARLY_ZERO, MAX_AGE

        data_list = self.filter_data(data_type)

        M,C = uninformative_prior_gp()

        age = []
        val = []
        V = []
        for d in data_list:
            if d['age_end'] == MISSING:
                d['age_end'] = MAX_AGE-1

            age.append(.5 * (d['age_start'] + d['age_end']))
            val.append(self.value_per_1(d) + .00001)
            V.append(self.se_per_1(d) ** 2.)

            if len(data_list) > 0:
                gp.observe(M, C, age, mc.logit(val), V)

        # use prior to set estimate near zero as requested
        near_zero = min(1., val)**2
        if near_zero == 1.:
            near_zero = 1e-9

        for prior_str in self.get_priors(data_type).split('\n'):
            prior = prior_str.split()
            if len(prior) > 0 and prior[0] == 'zero':
                age_start = int(prior[1])
                age_end = int(prior[2])

                gp.observe(M, C, range(age_start, age_end+1), mc.logit(near_zero), [0.])

        x = self.get_estimate_age_mesh()
        normal_approx_vals = mc.invlogit(M(x))

        self.set_initial_value(data_type, normal_approx_vals)

    
def get_disease_model(disease_model_id):
    """
    fetch specificed disease model data from
    dismod server given in settings.py
    """
    
    twc.go(DISMOD_LOGIN_URL)
    twc.fv('1', 'username', DISMOD_USERNAME)
    twc.fv('1', 'password', DISMOD_PASSWORD)
    twc.submit()
    twc.url('accounts/profile')

    twc.go(DISMOD_DOWNLOAD_URL % disease_model_id)

    return DiseaseJson(twc.show())

def post_disease_model(disease):
    """
    fetch specificed disease model data from
    dismod server given in settings.py
    """
    
    twc.go(DISMOD_UPLOAD_URL)
    twc.fv('1', 'model_json', disease.to_json())
    twc.submit()

    return twc.browser.get_url()
