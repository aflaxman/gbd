import numpy as np
import twill.commands as twc
import simplejson as json
from dismod3.settings import *

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
        self.set_key_by_type('priors', type)

    def get_estimate_age_mesh(self):
        return self.params.get('estimate_age_mesh', [])
    def set_estimate_age_mesh(self, mesh):
        self.params['estimate_age_mesh'] = list(mesh)

    def get_param_age_mesh(self):
        return self.params['param_age_mesh']
    def set_param_age_mesh(self, mesh):
        self.params['param_age_mesh'] = list(mesh)
    
    def filter_data(self, data_type):
        return [d for d in self.data if d['data_type'] == data_type]
    
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
