""" Generate a posterior estimate for a specific region, sex, and year
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

import numpy as np
import pymc as mc

import dismod3

def fit_posterior(id, region, sex, year):
    """ Fit posterior of specified region/sex/year for specified model

    Parameters
    ----------
    id : int
      The model id number for the job to fit
    region : str
      From dismod3.settings.gbd_regions, but clean()-ed
    sex : str, from dismod3.settings.gbd_sexes
    year : str, from dismod3.settings.gbd_years

    Example
    -------
    >>> import fit_posterior
    >>> fit_posterior.fit_posterior(2552, 'asia_east', 'male', '2005')
    """

    ## load model
    dm = dismod3.load_disease_model(id)


    ## separate out prevalence and relative-risk data
    prev_data = [d for d in dm.data if dismod3.gbd_disease_model.relevant_to(d, 'prevalence', region, year, sex)]
    rr_data = [d for d in dm.data if dismod3.gbd_disease_model.relevant_to(d, 'relative-risk', region, year, sex)]
    dm.data = [d for d in dm.data if not d in prev_data and not d in rr_data]


    ### setup the generic disease model (without prevalence data)
    import dismod3.gbd_disease_model as model
    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])
    dm.calc_effective_sample_size(dm.data)
    dm.vars = model.setup(dm, keys)


    ## override the birth prevalence prior, based on the withheld prevalence data
    logit_C_0 = dm.vars[dismod3.utils.gbd_key_for('bins', region, year, sex)]['initial'][2]
    assert len(prev_data) == 1, 'should be a single prevalance datum'
    d = prev_data[0]

    mu_logit_C_0 = mc.logit(dm.value_per_1(d)+dismod3.NEARLY_ZERO)
    lb, ub = dm.bounds_per_1(d)
    sigma_logit_C_0 = (mc.logit(ub+dismod3.NEARLY_ZERO) - mc.logit(lb+dismod3.NEARLY_ZERO)) / (2 * 1.96)
    print 'mu_C_0_pri:', mc.invlogit(mu_logit_C_0)
    print 'ui_C_0_pri:', lb, ub

    # override the excess-mortality, based on the relative-risk data
    if len(rr_data) > 0:
        d = rr_data[0]
        mu_rr = dm.value_per_1(d)
        sigma_rr = dm.se_per_1(d)

        log_f = dm.vars[dismod3.utils.gbd_key_for('excess-mortality', region, year, sex)]['age_coeffs']
        m_all = dm.vars[dismod3.utils.gbd_key_for('all-cause_mortality', region, year, sex)]
        mu_log_f = np.log((mu_rr-1) * m_all)
        sigma_log_f = 1 / ((mu_rr-1) * m_all) * np.log(sigma_rr * m_all)
    
    ### fit the model
    ## first generate decent initial conditions
    #model.fit(dm, method='map', keys=keys, verbose=dismod3.settings.ON_SGE)

    ## then sample the posterior via MCMC
    dm.mcmc = mc.MCMC(dm.vars)
    dm.mcmc.use_step_method(SampleFromNormal, logit_C_0, mu=mu_logit_C_0, tau=sigma_logit_C_0**-2)
    if len(rr_data) > 0:
        param_mesh = log_f.parents['param_mesh']
        log_f_mesh = log_f.parents['gamma_mesh']
        dm.mcmc.use_step_method(SampleFromNormal, log_f_mesh, mu=mu_log_f[param_mesh], tau=sigma_log_f[param_mesh]**-2)
    dm.mcmc.sample(1000, verbose=dismod3.settings.ON_SGE)

    print 'mu_C_0_post:', mc.invlogit(logit_C_0.stats()['mean'])
    print 'ui_C_0_post:', mc.invlogit(logit_C_0.stats()['95% HPD interval'])

    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    model.save_fit(dm, keys)
    dm.save('dm-%d-posterior-%s.json' % (id, dismod3.utils.gbd_key_for('all', region, sex, year)), keys_to_save=keys)

    return dm

class SampleFromNormal(mc.Gibbs):
    def __init__(self, stochastic, mu, tau, proposal_sd=None, verbose=None):
        mc.Gibbs.__init__(self, stochastic, verbose=verbose)
        self.mu = mu
        self.tau = tau
        
    def step(self):
        self.stochastic.value = mc.rnormal(self.mu, self.tau)

def main():
    import optparse

    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-s', '--sex', default='male',
                      help='only estimate given sex (valid settings ``male``, ``female``, ``all``)')
    parser.add_option('-y', '--year', default='2005',
                      help='only estimate given year (valid settings ``1990``, ``2005``)')
    parser.add_option('-r', '--region', default='australasia',
                      help='only estimate given GBD Region')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    dm = fit_posterior(id, options.region, options.sex, options.year)
    return dm

if __name__ == '__main__':
    dm = main()

