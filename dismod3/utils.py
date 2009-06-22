import pylab as pl
import numpy as np
import pymc as mc
from pymc import gp

from dismod3.settings import *

def trim(x, a, b):
    return np.maximum(a, np.minimum(b, x))

def const_func(x, c):
    """
    useful function for defining a non-informative
    prior on a Gaussian process
    >>> const_func([1,2,3], 17.0)
    [17., 17., 17.]
    """
    return np.zeros(len(x)) + c

def uninformative_prior_gp(c=-10.,  diff_degree=2., amp=100., scale=200.):
    """
    return mean and covariance objects for an uninformative prior on
    the age-specific rate
    """
    M = gp.Mean(const_func, c=c)
    C = gp.Covariance(gp.matern.euclidean, diff_degree=diff_degree,
                      amp=amp, scale=scale)

    return M, C

def spline_interpolate(in_mesh, values, out_mesh):
    from scipy.interpolate import interp1d
    f = interp1d(in_mesh, values, kind='linear')
    return f(out_mesh)

# def gp_interpolate(in_mesh, values, out_mesh):
#     """
#     interpolate a set of values given at
#     points on in_mesh to find values on
#     out_mesh.
#     """
#     M, C = uninformative_prior_gp()
#     gp.observe(M, C, in_mesh, values)
#     return M(out_mesh)

def interpolate(in_mesh, values, out_mesh):
    """
    wrapper so that it is only necessary to
    make one change to try different interpolation
    methods
    """
    return spline_interpolate(in_mesh, values, out_mesh)

def rate_for_range(raw_rate,age_indices,age_weights):
    """
    calculate rate for a given age-range,
    using the age-specific population numbers
    given by entries in years t0-t1 of pop_table,
    for given country and sex

    age_indices is a list of which indices of the raw rate
    should be used in the age weighted average (pre-computed
    because this is called in the inner loop of the mcmc)
    """
    age_adjusted_rate = np.sum(raw_rate[age_indices]*age_weights)/np.sum(age_weights)
    return age_adjusted_rate

def gbd_keys(type_list=stoch_var_types,
             region_list=gbd_regions,
             year_list=gbd_years,
             sex_list=gbd_sexes):
    """ Make a list of gbd keys for the type, region, year, and sex
    specified

    Parameters
    ----------
    type_list : list, optional, subset of ['incidence', 'remission', 'case-fatality']
    region_list : list, optional, subset of 21 GBD regions
    year_list : list, optional, subset of ['1990', '2005']
    sex_list : list, optional, subset of ['male', 'female']

    Results
    -------
    A list of gbd keys corresponding to all combinations of list
    items.
    """
    key_list = []

    # special case: prevalence is controlled by incidence, remission,
    # and case-fatality
    if type_list == [ 'prevalence' ]:
        types = [clean(t) for t in output_data_types]
        
    for t in type_list:
        for r in region_list:
            for y in year_list:
                for s in sex_list:
                    key_list.append(gbd_key_for(t, r, y, s))
    return key_list

def clean(str):
    """ Return a 'clean' version of a string, suitable for using as a hash
    string or a class attribute.
    """
    str = str.strip()
    str = str.lower()
    str = str.replace(',', '')
    str = str.replace('/', '_')
    str = str.replace(' ', '_')
    str = str.replace('(', '')
    str = str.replace(')', '')
    return str

def gbd_key_for(type, region, year, sex):
    """ Make a human-readable string that can be used as a key for
    storing estimates for the given type/region/year/sex.
    """
    return KEY_DELIM_CHAR.join([clean(type), clean(region),
                                str(year), clean(sex)])
    
def type_region_year_sex_from_key(key):
    return key.split(KEY_DELIM_CHAR)

def indices_for_range(age_mesh, age_start, age_end):
    return [ ii for ii, a in enumerate(age_mesh) if a >= age_start and a <= age_end ]

def prior_vals(dm, type):
    """ Estimate the prior distribution on param_age_mesh for a particular type

    Parameters
    ----------
    dm : DiseaseJson
    type : str, one of 'prevalence', 'incidence', 'remission', 'case-fatality'

    Results
    -------
    ages, vals : a list of ages and a list of estimated values

    TODO
    ----
    Instead of setting up the stochastics here, use the ones from the
    beta_binomial_model.setup (DRY)
    """
    param_mesh = dm.get_param_age_mesh()
    est_mesh = dm.get_estimate_age_mesh()

    logit_vals = mc.Normal('logit_vals',
                           mu=-5 * np.ones(len(param_mesh)),
                           tau=1.e-2,
                           value=-1 * np.ones(len(param_mesh)))

    conf = mc.Normal('conf', mu=100.0, tau=1./(30.)**2)

    @mc.deterministic
    def vals(logit_vals=logit_vals):
        return interpolate(param_mesh, mc.invlogit(logit_vals), est_mesh)

    prior_str = dm.get_global_priors(type)
    priors = generate_prior_potentials(prior_str, est_mesh, vals, conf)
    m = mc.MCMC([logit_vals, vals, priors])
    m.use_step_method(mc.AdaptiveMetropolis, logit_vals)
    m.sample(1000)

    m = mc.MCMC([conf])
    m.sample(100,50)

    return est_mesh, vals.value, conf.stats()['mean']

def prior_dict_to_str(pd):
    """ Generate a string suitable for passing to generate_prior_potentials
    from a prior dictionary

    Input
    -----
    pd : dict

    Notes
    -----
    This is a bit brittle, and a lot of duplicated code.  It should be rethought one day.
    """
    prior_str = ','

    smooth_str = {
        'No Prior': '',
        'Slightly': 'smooth 1,',
        'Moderately': 'smooth 10,',
        'Very': 'smooth 100,',
        }

    conf_str = {
        'No Prior': '',
        'Slightly': 'confidence 10 1,',
        'Moderately': 'confidence 50 1,',
        'Very': 'confidence 100 1,',
        }

    prior_str += smooth_str[pd.get('smoothness', 'No Prior')]
    prior_str += conf_str[pd.get('confidence', 'No Prior')]

    v = int(pd.get('zero_range', {}).get('age_before',0)) - 1
    if v >= 0:
        prior_str += 'zero 0 %d, ' % v
    
    v = int(pd.get('zero_range', {}).get('age_after',100)) + 1
    if v <= 100:
        prior_str += 'zero %d 100, ' % v

    v = float(pd.get('peak_bounds', {}).get('upper',1.0))
    if v < 1.:
        prior_str += 'max_at_most %f, ' % v

    v = float(pd.get('peak_bounds', {}).get('lower',0.0))
    if v > 0.:
        prior_str += 'max_at_least %f, ' % v
    
    return prior_str

def generate_prior_potentials(prior_str, age_mesh, rate, confidence_stoch=None):
    """
    return a list of potentials that model priors on the rate_stoch

    prior_str may have entries in the following format:
      smooth <tau> [<age_start> <age_end>]
      zero <age_start> <age_end>
      confidence <mean> <tau>
      increasing <age_start> <age_end>
      decreasing <age_start> <age_end>
      convex_up <age_start> <age_end>
      convex_down <age_start> <age_end>
      unimodal <age_start> <age_end>
      value <mean> <tau> [<age_start> <age_end>]
      max_at_least <value>
      max_at_most <value>
            
    for example: 'smooth .1, zero 0 5, zero 95 100'

    age_mesh[i] indicates what age the value of rate[i] corresponds to

    confidence_stoch can be an additional stochastic variable that is used
    in the beta-binomial model, but it is not required
    """

    def derivative_sign_prior(rate, prior, deriv, sign):
        age_start = int(prior[1])
        age_end = int(prior[2])
        age_indices = indices_for_range(age_mesh, age_start, age_end)
        @mc.potential(name='deriv_sign_{%d,%d,%d,%d}^%s' % (deriv, sign, age_start, age_end, rate))
        def deriv_sign_rate(f=rate,
                            age_indices=age_indices,
                            tau=1.e7,
                            deriv=deriv, sign=sign):
            df = np.diff(f[age_indices], deriv)
            return -tau * np.dot(df**2, (sign * df < 0))
        return [deriv_sign_rate]

    priors = []
    for line in prior_str.split(PRIOR_SEP_STR):
        prior = line.strip().split()
        if len(prior) == 0:
            continue
        if prior[0] == 'smooth':
            tau_smooth_rate = float(prior[1])

            if len(prior) == 4:
                age_start = int(prior[2])
                age_end = int(prior[3])
            else:
                age_start = 0
                age_end = MAX_AGE
            age_indices = indices_for_range(age_mesh, age_start, age_end)
                
            @mc.potential(name='smooth_{%d,%d}^%s' % (age_start, age_end, rate))
            def smooth_rate(f=rate, age_indices=age_indices, tau=tau_smooth_rate):
                return mc.normal_like(np.diff(np.log(np.maximum(NEARLY_ZERO, f[age_indices]))), 0.0, tau)
            priors += [smooth_rate]

        elif prior[0] == 'zero':
            tau_zero_rate = 1./(1e-4)**2
            
            age_start = int(prior[1])
            age_end = int(prior[2])
            age_indices = indices_for_range(age_mesh, age_start, age_end)
                               
            @mc.potential(name='zero_{%d,%d}^%s' % (age_start, age_end, rate))
            def zero_rate(f=rate, age_indices=age_indices, tau=tau_zero_rate):
                return mc.normal_like(f[age_indices], 0.0, tau)
            priors += [zero_rate]

        elif prior[0] == 'value':
            val = float(prior[1])
            tau = float(prior[2])

            if len(prior) == 4:
                age_start = int(prior[3])
                age_end = int(prior[4])
            else:
                age_start = 0
                age_end = MAX_AGE
            age_indices = indices_for_range(age_mesh, age_start, age_end)
                               
            @mc.potential(name='value_{%2f,%2f,%d,%d}^%s' \
                              % (val, tau, age_start, age_end, rate))
            def val_for_rate(f=rate, age_indices=age_indices, val=val, tau=tau):
                return mc.normal_like(f[age_indices], val, tau)
            priors += [val_for_rate]

        elif prior[0] == 'confidence':
            # prior only affects beta_binomial_rate model
            continue # ignore for now, while testing
            if not confidence_stoch:
                continue

            mu = float(prior[1])
            tau = float(prior[2])

            @mc.potential(name='prior_%s' % confidence_stoch)
            def confidence_potential(f=confidence_stoch, mu=mu, tau=tau):
                return mc.normal_like(f, mu, tau)
            priors += [confidence_potential]

        elif prior[0] == 'increasing':
            priors += derivative_sign_prior(rate, prior, deriv=1, sign=1)
        elif prior[0] == 'decreasing':
            priors += derivative_sign_prior(rate, prior, deriv=1, sign=-1)
        elif prior[0] == 'convex_down':
            priors += derivative_sign_prior(rate, prior, deriv=2, sign=-1)
        elif prior[0] == 'convex_up':
            priors += derivative_sign_prior(rate, prior, deriv=2, sign=1)

        elif prior[0] == 'unimodal':
            age_start = int(prior[1])
            age_end = int(prior[2])
            age_indices = indices_for_range(age_mesh, age_start, age_end)

            @mc.potential(name='unimodal_{%d,%d}^%s' % (age_start, age_end, rate))
            def unimodal_rate(f=rate, age_indices=age_indices, tau=1000.):
                df = np.diff(f[age_indices])
                sign_changes = pl.find((df[:-1] > NEARLY_ZERO) & (df[1:] < -NEARLY_ZERO))
                sign = np.ones(len(age_indices)-2)
                if len(sign_changes) > 0:
                    change_age = sign_changes[len(sign_changes)/2]
                    sign[change_age:] = -1.
                return -tau*np.dot(np.abs(df[:-1]), (sign * df[:-1] < 0))
            priors += [unimodal_rate]

        elif prior[0] == 'max_at_least':
            val = float(prior[1])

            @mc.potential(name='max_at_least{%f}^%s' % (val, rate))
            def max_at_least(f=rate, at_least=val, tau=1000./val**2):
                cur_max = np.max(f)
                return -tau * (cur_max - at_least)**2 * (cur_max < at_least)
            priors += [max_at_least]

        elif prior[0] == 'max_at_most':
            val = float(prior[1])

            @mc.potential(name='max_at_most{%f}^%s' % (val, rate))
            def max_at_most(f=rate, at_most=val, tau=1000./val**2):
                cur_max = np.max(f)
                return -tau * (cur_max - at_most)**2 * (cur_max > at_most)
            priors += [max_at_most]

        else:
            raise KeyError, 'Unrecognized prior: %s' % prior_str

    return priors






