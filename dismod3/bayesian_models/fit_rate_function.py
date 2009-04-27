import inspect

import numpy as np
import pymc as mc
from pymc import gp
from probabilistic_utils import uninformative_prior_gp, NEARLY_ZERO, MAX_AGE

MISSING = -99

import probabilistic_utils
from probabilistic_utils import trim

import beta_binomial_rate as rate_model
#import urbanicity_covariate_rate as rate_model

MAP_PARAMS = {
    'bfgs': [ 500, 'fmin_l_bfgs_b'],
    'powell': [ 500, 'fmin_powell'],
    'most accurate': [ 1500, 'fmin_powell'],
    'fast': [ 50, 'fmin_powell'],
    'testing fast': [ 1, 'fmin' ],
    }

MCMC_PARAMS = {
    'over accurate': [500, 200, 10000],
    'most accurate': [100, 50, 5000],
    'fast': [500, 10, 5000],
    'testing fast': [5, 5, 5],
    }

def initialize_model(asrf):
    # store the rate model code in the asrf for future reference
    asrf.fit['rate_model'] = inspect.getsource(rate_model)
    asrf.fit['out_age_mesh'] = range(probabilistic_utils.MAX_AGE)

    # do normal approximation first, to generate a good starting point
    M,C = probabilistic_utils.normal_approx(asrf)

    # define the model variables
    rate_model.setup_rate_model(asrf)

def map_fit(asrf, speed='most accurate'):
    """
    The Maximum A Posteriori (MAP) fit of the model is a point
    estimate of the model parameters, which is found using numerical
    optimization to attempt to maximize the posterior liklihood.
    Since this is a local optimization method, it might not find the
    global optimum.
    """
    initialize_model(asrf)
    
    map = mc.MAP(asrf.vars)
    iterlim, method = MAP_PARAMS[speed]
    print "searching for maximum likelihood point estimate (%s)" % method
    map.fit(verbose=1, iterlim=iterlim, method=method)

    probabilistic_utils.save_map(asrf)

def mcmc_fit(asrf, speed='most accurate'):
    """
    The Markov Chain Monte Carlo (MCMC) fit of the model works by
    making successive draws of the model parameters from the posterior
    distribution.  This provides confidence intervals, and should be
    more robust against local maxima in the posterior liklihood.  But
    the question is, did the chain run for long enough to mix?
    """
    map_fit(asrf, speed)
    
    print "drawing samples from posterior distribution (MCMC) (speed: %s)" % speed
    mcmc = mc.MCMC(asrf.vars)
    if asrf.vars.has_key('beta_binom_stochs'):
        mcmc.use_step_method(mc.AdaptiveMetropolis, asrf.vars['beta_binom_stochs'], verbose=0)
    if asrf.vars.has_key('logit(Erf_%d)'%asrf.id):
        mcmc.use_step_method(mc.AdaptiveMetropolis, asrf.vars['logit(Erf_%d)'%asrf.id], verbose=0)

    trace_len, thin, burn = MCMC_PARAMS[speed]
    mcmc.sample(trace_len*thin+burn, burn, thin, verbose=1)

    probabilistic_utils.save_mcmc(asrf)


    
import twill.commands as twc
import simplejson as json
from dismod3.settings import *


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

    return json.loads(twc.show())

def post_disease_model(disease_model):
    """
    fetch specificed disease model data from
    dismod server given in settings.py
    """
    
    twc.go(DISMOD_UPLOAD_URL)
    twc.fv('1', 'model_json', json.dumps(disease_model))
    twc.submit()

    return twc.browser.get_url()


def fit(disease_model, data_type='prevalence data'):
    """
    fit a single data_type from the model
    """
    dm = get_disease_model(disease_model)

    # filter out all data with type != data_type
    dm['data'] = [d for d in dm['data'] if d['data_type'] == data_type]

    # store the probabilistic model code for future reference
    dm['params'].update(bayesian_model=inspect.getsource(rate_model),
                     age_mesh=[0.0, 0.5, 3.0, 10.0, 20.0, 30.0, 40.0,
                               50.0, 60.0, 70.0, 80.0, 90.0, 100.0],
                     out_age_mesh=range(probabilistic_utils.MAX_AGE),
                     map={})

    # do normal approximation first, to generate a good starting point
    fit_normal_approx(dm, data_type)
    
    # define the model variables
    bayesian_model = setup_rate_model(dm, data_type)

    map = mc.MAP(bayesian_model)
    map.fit(method='fmin_powell', iterlim=500, tol=.001, verbose=1)
    dm['params']['map'][data_type] = list(bayesian_model.rate.value)
    
    post_disease_model(dm)
    return dm, bayesian_model


def setup_rate_model(dm, data_type):
    #############################################################################
    # set up the age-specific Beta stochastic variables
    #
    MIN_CONFIDENCE = 1
    MAX_CONFIDENCE = 100000
    
    initial_value = np.array(dm['params']['normal_approx'][data_type])
    mesh = dm['params']['age_mesh']
    out_mesh = dm['params']['out_age_mesh']
    rate_str = data_type.replace('data','')

    logit_rate = mc.Normal('logit(%s)' % rate_str,
                           mu=np.zeros(len(mesh)),
                           tau=1.e-2,
                           value=mc.logit(initial_value[mesh]),
                           verbose=0)
    @mc.deterministic(name=rate_str)
    def rate(logit_rate=logit_rate):
        return probabilistic_utils.interpolate(mesh, mc.invlogit(logit_rate), out_mesh)
    
    confidence = mc.Normal('conf_%s' % rate_str, mu=1000.0, tau=1./(300.)**2)
    
    @mc.deterministic(name='alpha_%s' % rate_str)
    def alpha(rate=rate, confidence=confidence):
        return rate * trim(confidence, MIN_CONFIDENCE, MAX_CONFIDENCE)

    @mc.deterministic(name='beta_%s' % rate_str)
    def beta(rate=rate, confidence=confidence):
        return (1. - rate) * trim(confidence, MIN_CONFIDENCE, MAX_CONFIDENCE)

    ########################################################################
    # set up stochs for the priors and observed data
    #
    prior_str = dm['params'].get('priors', {}).get(data_type, '')
    priors = generate_prior_potentials(prior_str, rate, confidence)

    logit_p_stochs = []
    p_stochs = []
    beta_potentials = []
    observed_rates = []
    for d in dm['data']:
        if d['data_type'] != data_type:
            continue
        id = d['id']
        
        # ensure all rate data is valid
        # TODO: raise exceptions to have users fix an errors
        d.update(value=trim(d['value'], NEARLY_ZERO, 1.-NEARLY_ZERO),
                 standard_error=max(d['standard_error'], 0.0001))

        logit_p = mc.Normal('logit(p_%d)' % id, 0., 1/(10.)**2,
                            value=mc.logit(d['value'] + NEARLY_ZERO),
                            verbose=0)

        p = mc.InvLogit('p_%d' % id, logit_p)

        @mc.potential(name='beta_potential_%d' % id)
        def potential_p(p=p,
                        alpha=alpha, beta=beta,
                        a0=d['age_start'], a1=d['age_end'],
                        age_weights=d['age_weights']):
            a = probabilistic_utils.rate_for_range(alpha, a0, a1, age_weights)
            b = probabilistic_utils.rate_for_range(beta, a0, a1, age_weights)
            return mc.beta_like(trim(p, NEARLY_ZERO, 1. - NEARLY_ZERO), a, b)

        denominator = max(100., d['value'] * (1 - d['value']) / d['standard_error']**2.)
        numerator = d['value'] * denominator
        obs = mc.Binomial("data_%d" % id, value=numerator, n=denominator, p=p, observed=True)

        logit_p_stochs.append(logit_p)
        p_stochs.append(p)
        beta_potentials.append(potential_p)
        observed_rates.append(obs)
        
    return mc.Model(locals())


def extract_units(d):
    """
    d is a data hash which might include
    the key 'units', which is a decription
    of the units for this datum.

    return the float that d['value'] should
    be multiplied to make the units per 1.0
    """
    return 1. / float(d['units'].split()[-1])

def fit_normal_approx(dm, data_type):
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
    param_hash = dm['params']
    data_list = [d for d in dm['data'] if d['data_type'] == data_type]

    M,C = uninformative_prior_gp()

    age = []
    val = []
    V = []
    for d in data_list:
        scale = extract_units(d)

        if d['age_end'] == MISSING:
            d['age_end'] = MAX_AGE-1

        d['standard_error'] /= scale
        if d['standard_error'] == 0.:
            d['standard_error'] = .001

        d['value'] *= scale
        d['units'] = 'per 1.0'

        age.append(.5 * (d['age_start'] + d['age_end']))
        val.append(d['value'] + .00001)
        V.append((d['standard_error']) ** 2.)

    if len(data_list) > 0:
        gp.observe(M, C, age, mc.logit(val), V)

    # use prior to set estimate near zero as requested
    near_zero = min(1., val)**2
    if near_zero == 1.:
        near_zero = 1e-9
        
    for prior_str in param_hash.get('priors', '').split('\n'):
        prior = prior_str.split()
        if len(prior) > 0 and prior[0] == 'zero':
            age_start = int(prior[1])
            age_end = int(prior[2])

            gp.observe(M, C, range(age_start, age_end+1), mc.logit(near_zero), [0.])
        
    x = param_hash['out_age_mesh']
    normal_approx_vals = mc.invlogit(M(x))

    if not param_hash.has_key('normal_approx'):
        param_hash['normal_approx'] = {}

    param_hash['normal_approx'][data_type] = list(normal_approx_vals)


def generate_prior_potentials(prior_str, rate, confidence):
    """
    return a list of potentials that model priors on the rate_stoch
    prior_str may have lines in the following format:
      smooth <tau> <age_start> <age_end>
      zero <age_start> <age_end>
      confidence <mean> <tau>
      increasing <age_start> <age_end>
      decreasing <age_start> <age_end>
      convex_up <age_start> <age_end>
      convex_down <age_start> <age_end>
      unimodal <age_start> <age_end>
    
    for example: 'smooth .1 \n zero 0 5 \n zero 95 100'
    """

    def derivative_sign_prior(rate, prior, deriv, sign):
        age_start = int(prior[1])
        age_end = int(prior[2])
        @mc.potential(name='deriv_sign_{%d,%d,%d,%d}^%s' % (deriv, sign, age_start, age_end, rate))
        def deriv_sign_rate(f=rate,
                            age_start=age_start, age_end=age_end,
                            tau=1000.,
                            deriv=deriv, sign=sign):
            df = np.diff(f[age_start:(age_end+1)], deriv)
            return -tau * np.dot(df**2, (sign * df < 0))
        return [deriv_sign_rate]

    priors = []
    
    for line in prior_str.split('\n'):
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
                
            @mc.potential(name='smooth_{%d,%d}^%s' % (age_start, age_end, rate))
            def smooth_rate(f=rate, age_start=age_start, age_end=age_end, tau=tau_smooth_rate):
                return mc.normal_like(np.diff(np.log(np.maximum(NEARLY_ZERO, f[range(age_start, age_end)]))), 0.0, tau)
            priors += [smooth_rate]

        elif prior[0] == 'zero':
            tau_zero_rate = 1./(1e-4)**2
            
            age_start = int(prior[1])
            age_end = int(prior[2])
                               
            @mc.potential(name='zero_{%d,%d}^%s' % (age_start, age_end, rate))
            def zero_rate(f=rate, age_start=age_start, age_end=age_end, tau=tau_zero_rate):
                return mc.normal_like(f[range(age_start, age_end+1)], 0.0, tau)
            priors += [zero_rate]

        elif prior[0] == 'confidence':
            # prior only affects beta_binomial_rate model
            if not confidence:
                continue

            mu = float(prior[1])
            tau = float(prior[2])

            @mc.potential(name='prior_%s' % confidence)
            def confidence(f=confidence, mu=mu, tau=tau):
                return mc.normal_like(f, mu, tau)
            priors += [confidence]

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
            @mc.potential(name='unimodal_{%d,%d}^%s' % (age_start, age_end, rate))
            def unimodal_rate(f=rate, age_start=age_start, age_end=age_end, tau=1000.):
                df = np.diff(f[age_start:(age_end + 1)])
                sign_changes = pl.find((df[:-1] > NEARLY_ZERO) & (df[1:] < -NEARLY_ZERO))
                sign = np.ones(age_end-age_start-1)
                if len(sign_changes) > 0:
                    change_age = sign_changes[len(sign_changes)/2]
                    sign[change_age:] = -1.
                return -tau*np.dot(np.abs(df[:-1]), (sign * df[:-1] < 0))
            priors += [unimodal_rate]

        else:
            raise KeyException, 'Unrecognized prior: %s' % prior_str

    return priors




import pylab as pl

def plot_disease_model(dm):
    
    # divide up disease_model data by data_type
    #import pdb; pdb.set_trace()
    data_by_type = {}
    for d in dm['data']:
        data_by_type[d['data_type']] = data_by_type.get(d['data_type'], []) + [d]

    types = set(data_by_type.keys()) | set(dm['params'].get('map', {}).keys())
    cnt = max(1, len(types))
    cols = min(2, cnt)
    rows = int(np.ceil(float(cnt) / float(cols)))

    subplot_width = 6
    subplot_height = 4
    
    clear_plot(width=subplot_width*cols,height=subplot_height*rows)
    for ii, type in enumerate(types):
        data = data_by_type.get(type, [])
        
        pl.subplot(rows, cols, ii + 1)
        plot_intervals(data, fontsize=12, alpha=.5)
        plot_normal_approx(dm, type)
        plot_map_fit(dm, type)
        #plot_mcmc_fit(dm, type)
        plot_truth(dm, type)
        #plot_prior(dm, type)
        label_plot('%s (id=%d)' % (type, dm['params']['id']), fontsize=10)
        
        max_rate = np.max([.0001] + [d['value']*extract_units(d) for d in data])
        xmin = 0.
        xmax = 85.
        ymin = 0.
        ymax = 1.25*max_rate
        pl.axis([xmin, xmax, ymin, ymax])

        if ii % cols != 0:
            pl.ylabel('')
        if (ii + cols) < cnt:
            pl.xlabel('')

def plot_intervals(data, alpha=.75, color=(.0,.5,.0), text_color=(.0,.3,.0), fontsize=12):
    """
    use matplotlib plotting functions to render transparent
    rectangles on the current figure representing each
    piece of Data
    """
    for d in data:
        if d['age_end'] == MISSING:
            d['age_end'] = MAX_AGE

        scale = extract_units(d)
        val = d['value']*scale
        lower_ci = max(0., val - 1.98 * d['standard_error'] * scale)
        upper_ci = min(1., val + 1.98 * d['standard_error'] * scale)
        pl.plot([.5 * (d['age_start']+d['age_end']+1)]*2,
                [lower_ci, upper_ci],
                color=color, alpha=alpha, linewidth=1)
        pl.plot(np.array([d['age_start'], d['age_end']+1.]),
                np.array([val, val]),
                color=color, alpha=alpha, linewidth=5,
                )
    
def plot_fit(dm, fit_name, data_type, **params):
    fit = dm['params'].get(fit_name, {}).get(data_type)
    if fit:
        pl.plot(dm['params'].get('out_age_mesh'), fit, **params)

def plot_normal_approx(dm, type):
    plot_fit(dm, 'normal_approx', type, color='blue', alpha=.5)

def plot_truth(dm, type):
    plot_fit(dm, 'truth', type, linestyle='dotted', color='green', alpha=.95, linewidth=2)

def plot_map_fit(dm, type, **params):
    default_params = {'color': 'blue',
                      #'linestyle': 'dashed',
                      'linewidth': 2,
                      'alpha': .9,
                      }
    default_params.update(**params)
    plot_fit(dm, 'map', type, **default_params)

def plot_mcmc_fit(rf, detailed_legend=False, color=(.2,.2,.2)):
    try:
        x = np.concatenate((rf.fit['out_age_mesh'], rf.fit['out_age_mesh'][::-1]))
        y = np.concatenate((rf.fit['mcmc_lower_cl'], rf.fit['mcmc_upper_cl'][::-1]))

        pl.fill(x, y, facecolor='.2', edgecolor=color, alpha=.5)

        mcmc_avg = rf.fit['mcmc_mean']

        if detailed_legend:
            label = str(rf.region)
            color = np.random.rand(3)
        else:
            label = 'MCMC Fit'
            color = color

        pl.plot(rf.fit['out_age_mesh'], mcmc_avg, color=color, linewidth=3, alpha=.75, label=label)
    except (KeyError, ValueError):
        pass
        #pl.figtext(0.4,0.4, 'No MCMC Fit Found')


def clear_plot(width=4*1.5, height=3*1.5):
    fig = pl.figure(figsize=(width,height))
    pl.clf()
    return fig

def label_plot(title, **params):
    pl.xlabel('Age (years)', **params)
    pl.ylabel('Rate (per 1.0)', **params)
    pl.title(str(title), **params)
