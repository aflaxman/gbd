import numpy as np
import pymc as mc
from pymc import gp

import twill.commands as twc
import simplejson as json

from dismod3.settings import *

from bayesian_models.probabilistic_utils import trim, uninformative_prior_gp, \
    NEARLY_ZERO, MAX_AGE
MISSING = -99


class DiseaseJson:
    def __init__(self, json_str):
        dm = json.loads(json_str)
        self.params = dm['params']
        self.data = dm['data']

    def to_json(self):
        return json.dumps({'params': self.params, 'data': self.data})

    def set_region(self, region):
        """ Set the region of the disease model"""
        self.params['region'] = region
    def get_region(self):
        """ Get the region of the disease model"""
        return self.params.get('region', '')

    def set_ymax(self, val):
        """ Set the maximum y scale for plotting the disease model"""
        self.params['ymax'] = val
    def get_ymax(self):
        """ Get the maximum y scale for plotting the disease model"""
        return self.params.get('ymax', .1)
        
    def set_key_by_type(self, key, type, value):
        if not self.params.has_key(key):
            self.params[key] = {}
        self.params[key][type] = value
    def get_key_by_type(self, key, type, default=None):
        return self.params.get(key, {}).get(type, default)

    def get_initial_value(self, type):
        """ Return the initial value for estimate of a particular
        type, default to NEARLY_ZERO"""
        default_val = NEARLY_ZERO * np.ones(len(self.get_estimate_age_mesh()))
        return np.array(
            self.get_key_by_type('initial_value', type, default=default_val)
            )
    def set_initial_value(self, type, val):
        self.set_key_by_type('initial_value', type, list(val))
    def has_initial_value(self, type):
        return self.params.get('initial_value', {}).has_key(type)

    def get_map(self, type):
        return np.array(self.get_key_by_type('map', type, default=[]))
    def set_map(self, type, val):
        self.set_key_by_type('map', type, list(val))

    def get_mcmc(self, est_type, data_type):
        return np.array(self.get_key_by_type('mcmc_%s' % est_type, data_type, default=[]))
    def set_mcmc(self, est_type, data_type, val):
        self.set_key_by_type('mcmc_%s' % est_type, data_type, list(val))

    def clear_fit(self):
        """ Clear all estimates, fits, and stochastic vars

        Results
        -------
        disease_model.clear_fit() removes all the results of a fit, so
        that the model is ready for fitting fresh, for example with
        different new priors

        Example
        -------
        >>> import dismod3
        >>> dm = dismod3.get_disease_model(1)
        >>> dm.clear_fit()
        """
        for k in self.params.keys():
            if k == 'map' or k.find('mcmc_') >= 0:
                self.params.pop(k)

        if hasattr(self, 'vars'):
            delattr(self, 'vars')
            
        if hasattr(self, 'map'):
            delattr(self, 'map')
            
        if hasattr(self, 'mcmc'):
            delattr(self, 'mcmc')

    def get_units(self, type):
        return self.get_key_by_type('units', type)
    def set_units(self, type, units):
        self.set_key_by_type('units', type, units)

    def get_priors(self, type):
        """ Return the priors for estimates of given type"""
        return self.get_key_by_type('priors', type) or ''
    def set_priors(self, type, priors):
        """ Set the prior for data of a given type

        Parameters
        ----------
        type : str
          The type of data to which these priors apply
        priors : str
          The priors, see the prior generating function for details

        Notes
        -----
        Any stochastic variables are deleted, since they do not
        include the new priors
        """
        self.set_key_by_type('priors', type, priors)
        self.clear_fit()

    def get_estimate_age_mesh(self):
        return self.params.get('estimate_age_mesh', range(MAX_AGE))
    def set_estimate_age_mesh(self, mesh):
        """ Set the age mesh for the estimated age functions

        Parameters
        ----------
        mesh : list
          The estimate mesh.  Estimates for the prevalence, incidence,
          etc will be estimated at each point on this mesh.  For the
          generic disease model, the distance between consecutive mesh
          points must be one.

        Notes
        -----
        Any stochastic variables are deleted, since they need to be
        regenerated to reflect the new age mesh
        """
        self.params['estimate_age_mesh'] = list(mesh)
        self.clear_fit()
        
    def get_param_age_mesh(self):
        return self.params.get('param_age_mesh', range(0, MAX_AGE, 10))
    def set_param_age_mesh(self, mesh):
        """ Set the age mesh for the age functions control parameters

        Parameters
        ----------
        mesh : list
          The param mesh, the values at these mesh points will be
          linearly interpolated to form the age-specific function

        Notes
        -----
        To save time, the stochastic models use a linear interpolation
        of the points on the param age mesh to represent the age
        function.
        
        Any stochastic variables are deleted, since they need to be
        regenerated to reflect the new age mesh
        """
        self.params['param_age_mesh'] = list(mesh)
        self.clear_fit()
        
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
            if unit_str.strip()[0:4].lower() == 'per ':
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
        if self.params.get('initial_value',{}).has_key('all-cause mortality'):
            return self.get_initial_value('all-cause mortality')
        
        mortality_data = self.filter_data('all-cause mortality data')
        if len(mortality_data) == 0:
            return np.zeros(len(self.get_estimate_age_mesh()))
        else:
            self.fit_initial_estimate('all-cause mortality', mortality_data)
            return self.get_initial_value('all-cause mortality')

    def value_per_1(self, data):
        if data['value'] == MISSING:
            return MISSING

        return data['value'] * self.extract_units(data)

    def se_per_1(self, data):
        # TODO: extract se from ci if missing
        if data['standard_error'] == MISSING:
            return MISSING

        se = data['standard_error']
        se *= self.extract_units(data)

        return se
        

    def fit_initial_estimate(self, est_name, data_list):
        """ Find an initial estimate of the age-specific data

        Parameters
        ----------
        est_name : str
          The name of the estimate, as used in
          self.set_initial_value(est_name) and
          self.get_inital_value(est_name).

        data_list : list of data dicts
          The data to use for creating the initial estimate.
          
        Results
        -------
        The estimated parameter values are stored using
        self.set_initial_value(est_name, values)

        Example
        -------
        >>> import dismod3
        >>> dm = dismod3.get_disease_model(1)
        >>> dm.find_initial_estimate('prevalence', dm.filter_data('prevalence data'))
        >>> dm.get_initial_value('prevalence')
        
        Notes
        -----
        This estimate for an age-specific dataset is formed by using
        each datum to produce an estimate of the function value at a
        single age, and then saying that the logit of the true rate
        function is a gaussian process and these single age estimates
        are observations of this gaussian process.
        
        This is less valid and less accurate than using MCMC or MAP,
        but it can be much faster.  It is used to generate an initial
        value for the maximum-liklihood estimate.
        """
        from bayesian_models import probabilistic_utils
        from bayesian_models.probabilistic_utils import \
            trim, uninformative_prior_gp, NEARLY_ZERO, MAX_AGE

        M,C = uninformative_prior_gp()

        age = []
        val = []
        V = []
        for d in data_list:
            v = self.value_per_1(d)
            if v == MISSING:
                continue
            val.append(v + .00001)

            if d['age_end'] == MISSING:
                d['age_end'] = MAX_AGE-1
            age.append(.5 * (d['age_start'] + d['age_end']))

            se = self.se_per_1(d)
            if se == MISSING:
                V.append(.1)
            else:
                V.append((se+.00001) ** 2.)

            if len(data_list) > 0:
                gp.observe(M, C, age, mc.logit(val), V)

        # use prior to set estimate near zero as requested
        near_zero = min(1., val)**2
        if near_zero == 1.:
            near_zero = 1e-9

        for prior_str in self.get_priors(est_name).split('\n'):
            prior = prior_str.split()
            if len(prior) > 0 and prior[0] == 'zero':
                age_start = int(prior[1])
                age_end = int(prior[2])

                gp.observe(M, C, range(age_start, age_end+1), mc.logit(near_zero), [0.])

        x = self.get_estimate_age_mesh()
        normal_approx_vals = mc.invlogit(M(x))

        self.set_initial_value(est_name, normal_approx_vals)

    
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
