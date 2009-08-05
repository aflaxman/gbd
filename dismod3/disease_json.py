import numpy as np
import pymc as mc
from pymc import gp

import twill.commands as twc
import simplejson as json

from dismod3.settings import *
from dismod3.utils import debug, clean, trim, uninformative_prior_gp, prior_dict_to_str, NEARLY_ZERO, MAX_AGE, MISSING

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
        return self.params.get('ymax', 1.)

    def set_notes(self, val):
        """ Set notes for the disease model"""
        self.params['notes'] = val
    def get_notes(self):
        """ Get notes for the disease model"""
        return self.params.get('notes', '')
        
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
    def has_map(self, type):
        return self.params.get('map', {}).has_key(type)

    def get_truth(self, type):
        return np.array(self.get_key_by_type('truth', type, default=[]))
    def set_truth(self, type, val):
        self.set_key_by_type('truth', type, list(val))
    def has_truth(self, type):
        return self.params.get('truth', {}).has_key(type)

    def get_mcmc(self, est_type, data_type):
        return np.array(self.get_key_by_type('mcmc_%s' % est_type, data_type, default=[]))
    def set_mcmc(self, est_type, data_type, val):
        self.set_key_by_type('mcmc_%s' % est_type, data_type, list(val))

    def get_population(self, region):
        return np.array(self.get_key_by_type('population', region, default=None))
    def set_population(self, region, val):
        self.set_key_by_type('population', region, list(val))

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
        """ Return the priors for estimates of given type

        If no specific priors are found for the given type, default to the
        value in get_global_priors(type)
        """
        prior_str = self.get_key_by_type('priors', type, default='')
        if not prior_str:
            prior_str = self.get_global_priors(type)
        return prior_str
    
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

    def get_global_priors(self, type):
        """ Return the global priors that best match the specified type

        Since the type might be a key with the form
        'incidence+sub-saharan_africa_east+1990+female', return the
        first global prior who's key is found as a substring of ``type``

        Build and cache the global_priors_dict from the
        global_priors_json, if necessary.
        """
        if not hasattr(self, 'global_priors'):
            raw_dict = json.loads(self.params.get('global_priors_json', '{}'))
            self.global_priors = {'prevalence': {},
                                  'incidence': {},
                                  'remission': {},
                                  'case_fatality': {}}

            # reverse the order of the first and second level of keys in the raw_dict
            # this will be more convenient later
            for k1 in ['confidence', 'smoothness', 'zero_range', 'peak_bounds']:
                if not raw_dict.has_key(k1):
                    continue
                for k2 in raw_dict[k1]:
                    self.global_priors[k2][k1] = raw_dict[k1][k2]

            # deal with the dash vs underscore
            self.global_priors['case-fatality'] = self.global_priors['case_fatality']
            
            for k in self.global_priors:
                self.global_priors[k]['prior_str'] = prior_dict_to_str(self.global_priors[k])
        for k in self.global_priors:
            if clean(type).find(clean(k)) != -1:
                return self.global_priors[k]['prior_str']

        return ''

    def extract_params_from_global_priors(self):
        """ The global priors hash contains information on the age mesh,
        max y value, and additional notes, which should be stored
        somewhere more convenient
        """
        gp_dict = json.loads(self.params.get('global_priors_json', '{}'))
        self.set_param_age_mesh([int(a) for a in gp_dict['parameter_age_mesh']])
        self.set_ymax(float(gp_dict['y_maximum']))
        self.set_notes(gp_dict['note'])

    def set_empirical_prior(self, type, prior_dict):
        """ The empirical prior hash contains model-specific data for
        keyed by model parameter types"""
        self.set_key_by_type('empirical_prior', type, prior_dict)
    def get_empirical_prior(self, type):
        """ The empirical prior is a model specific dictionary"""
        return self.get_key_by_type('empirical_prior', type, default={})

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

        TODO: migrate to using 'radix', a number with no 'per '
        business

        This is hacky, so examples are best for now

        Example
        -------
        >>> dm.extract_units({})
        1.
        >>> dm.extract_units({'units': 'per 10'})
        .1
        >>> dm.extract_units({'units': '10'})
        .1
        >>> dm.extract_units({'units': 'bananas'})
        1.
        
        """
        try:
            unit_str = d.get('units', '1')
            unit_str = unit_str.replace('per ', '')
            unit_str = unit_str.replace(',', '')
            units = 1. / float(unit_str)
            return units
        except ValueError:
            debug('could not parse unit str: %s' % unit_str)
            return 1.


    def mortality(self, key='all-cause_mortality', data=None):
        """ Calculate the all-cause mortality rate for the
        region and sex of disease_model, and return it
        in an array corresponding to age_mesh

        Parameters
        ----------
        key : str, optional
          of the form 'all-cause_mortality+gbd_region+year+sex'
        data: list, optional
          the data list to extract all-cause mortality from
        """
        if self.params.get('initial_value',{}).has_key(key):
            return self.get_initial_value(key)

        if not data:
            data = self.filter_data('all-cause_mortality data')
        
        if len(data) == 0:
            return NEARLY_ZERO * np.ones(len(self.get_estimate_age_mesh()))
        else:
            self.fit_initial_estimate(key, data)
            return self.get_initial_value(key)

    def value_per_1(self, data):
        if data['value'] == MISSING:
            return MISSING

        return data['value'] * self.extract_units(data)

    def se_per_1(self, data):
        # TODO: extract se from ci if missing
        #   order of preference: use se if given, to calculate N s.t. stdev(Bi(N,p)) = se
        #   if no se, but upper_ui, lower_ui,
        #     if symmetric around p, find N assuming p is normal
        #     if non-symmetric, find N assuming logit(p) is normal
        if data['standard_error'] == MISSING:
            return MISSING

        se = data['standard_error']
        se *= self.extract_units(data)

        return se
        

    def fit_initial_estimate(self, key, data_list):
        """ Find an initial estimate of the age-specific data

        Parameters
        ----------
        key : str
          The name of the estimate, as used in
          self.set_initial_value(key) and
          self.get_inital_value(key).

        data_list : list of data dicts
          The data to use for creating the initial estimate.
          
        Results
        -------
        The estimated parameter values are stored using
        self.set_initial_value(key, values)

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
        from dismod3.logit_gp_step import LogitGPStep
        lr = mc.Normal('lr', -5 * np.ones(100), .1e-2)
        sm = LogitGPStep(lr, dm=self, key=key, data_list=data_list)
        x = self.get_estimate_age_mesh()
        normal_approx_vals = mc.invlogit(sm.M(x))
        self.set_initial_value(key, normal_approx_vals)
    
def get_disease_model(disease_model_id):
    """
    fetch specificed disease model data from
    dismod server given in settings.py
    """
    dismod_server_login()
    
    twc.go(DISMOD_DOWNLOAD_URL % disease_model_id)
    return DiseaseJson(twc.show())

def post_disease_model(disease):
    """
    fetch specificed disease model data from
    dismod server given in settings.py
    """
    dismod_server_login()

    # don't upload the disease data, since it is already on the server
    data = disease.data
    disease.data = []
    d_json = disease.to_json()
    disease.data = data

    twc.go(DISMOD_UPLOAD_URL)
    twc.fv('1', 'model_json', d_json)
    twc.submit()


    return twc.browser.get_url()

def get_job_queue():
    """
    fetch list of disease model jobs waiting to run from dismod server
    given in settings.py.
    """
    dismod_server_login()
    twc.go(DISMOD_LIST_JOBS_URL)
    return json.loads(twc.show())

def remove_from_job_queue(id):
    """
    remove a disease model from the job queue on the dismod server
    given in dismod3/settings.py
    """
    dismod_server_login()
    
    twc.go(DISMOD_REMOVE_JOB_URL)
    twc.fv('1', 'id', str(id))
    twc.submit()
    return twc.browser.get_url()
    

def dismod_server_login():
    """ login to the dismod server given in dismod3/settings.py."""
    
    twc.go(DISMOD_LOGIN_URL)
    twc.fv('1', 'username', DISMOD_USERNAME)
    twc.fv('1', 'password', DISMOD_PASSWORD)
    twc.submit()
    twc.url('accounts/profile')
    

