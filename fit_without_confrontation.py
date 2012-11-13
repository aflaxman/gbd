""" Generate a posterior estimate for a specific region, sex, and year
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

import numpy as np
import pymc as mc

import dismod3

def fit_without_confrontation(id, region, sex, year):
    """ Fit posterior of specified region/sex/year for specified model
    without trying to integrate conflicting sources of data

    Parameters
    ----------
    id : int
      The model id number for the job to fit
    region : str
      From dismod3.settings.gbd_regions, but clean()-ed
    sex : str, from dismod3.settings.gbd_sexes
    year : str, from dismod3.settings.gbd_years
    """

    ## load model
    dm = dismod3.load_disease_model(id)


    ## separate out prevalence and relative-risk data
    prev_data = [d for d in dm.data if dm.relevant_to(d, 'prevalence', region, year, sex)]
    rr_data = [d for d in dm.data if dm.relevant_to(d, 'relative-risk', region, year, sex)]
    dm.data = [d for d in dm.data if not d in prev_data and not d in rr_data]


    ### setup the generic disease model (without prevalence data)
    import dismod3.gbd_disease_model as model
    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])
    dm.calc_effective_sample_size(dm.data)
    dm.vars = model.setup(dm, keys)


    ## override the birth prevalence prior, based on the withheld prevalence data
    logit_C_0 = dm.vars[dismod3.utils.gbd_key_for('bins', region, year, sex)]['initial']['logit_C_0']
    assert len(prev_data) == 1, 'should be a single prevalance datum'
    d = prev_data[0]

    mu_logit_C_0 = mc.logit(dm.value_per_1(d)+dismod3.settings.NEARLY_ZERO)
    lb, ub = dm.bounds_per_1(d)
    sigma_logit_C_0 = (mc.logit(ub+dismod3.settings.NEARLY_ZERO) - mc.logit(lb+dismod3.settings.NEARLY_ZERO)) / (2 * 1.96)
    print 'mu_C_0_pri:', mc.invlogit(mu_logit_C_0)
    print 'ui_C_0_pri:', lb, ub

    # override the excess-mortality, based on the relative-risk data
    mu_rr = 1.01*np.ones(dismod3.settings.MAX_AGE)
    sigma_rr = .01*np.ones(dismod3.settings.MAX_AGE)
    for d in rr_data:
        mu_rr[d['age_start']:(d['age_end']+1)] = dm.value_per_1(d)
        sigma_rr[d['age_start']:(d['age_end']+1)] = dm.se_per_1(d)
    print 'mu_rr:', mu_rr.round(2)
    #print 'sigma_rr:', sigma_rr.round(2)

    log_f = dm.vars[dismod3.utils.gbd_key_for('excess-mortality', region, year, sex)]['age_coeffs']
    log_f_mesh = log_f.parents['gamma_mesh']
    param_mesh = log_f.parents['param_mesh']
    
    m_all = dm.vars[dismod3.utils.gbd_key_for('all-cause_mortality', region, year, sex)]
    mu_log_f = np.log((mu_rr-1) * m_all)
    sigma_log_f = 1 / ((mu_rr-1) * m_all) * sigma_rr * m_all
    print 'mu_log_f:', mu_log_f.round(2)[param_mesh]
    print 'sigma_log_f:', sigma_log_f.round(2)[param_mesh]
    
    ### fit the model using Monte Carlo simulation (shoehorned into the MCMC framework of PyMC)
    dm.mcmc = mc.MCMC(dm.vars)
    dm.mcmc.use_step_method(SampleFromNormal, logit_C_0, mu=mu_logit_C_0, tau=sigma_logit_C_0**-2)
    dm.mcmc.use_step_method(SampleFromNormal, log_f_mesh, mu=mu_log_f[param_mesh], tau=sigma_log_f[param_mesh]**-2)
    for stoch in dm.mcmc.stochastics:
        dm.mcmc.use_step_method(mc.NoStepper, stoch)
    dm.mcmc.sample(1000, verbose=dismod3.settings.ON_SGE)

    #print 'mu_C_0_post:', mc.invlogit(logit_C_0.stats()['mean']).round(2)
    #print 'ui_C_0_post:', mc.invlogit(logit_C_0.stats()['95% HPD interval']).round(2)
    #print 'mu_rr_post:', dm.vars[dismod3.utils.gbd_key_for('relative-risk', region, year, sex)]['rate_stoch'].stats()['mean'].round(2)
    print 'mu_log_f_mesh_post:', log_f_mesh.stats()['mean'].round(2)
    print 'mu_f_post:', dm.vars[dismod3.utils.gbd_key_for('excess-mortality', region, year, sex)]['rate_stoch'].stats()['mean'].round(2)


    for k in keys:
        t,r,y,s = dismod3.utils.type_region_year_sex_from_key(k)

        if t in ['incidence', 'prevalence', 'remission', 'excess-mortality', 'mortality', 'prevalence_x_excess-mortality']:
            dismod3.neg_binom_model.store_mcmc_fit(dm, k, dm.vars[k])

        elif t in ['relative-risk', 'duration', 'incidence_x_duration']:
            dismod3.normal_model.store_mcmc_fit(dm, k, dm.vars[k])

    from fit_posterior import save_country_level_posterior
    if str(year) == '2005':  # also generate 2010 estimates
        save_country_level_posterior(dm, region, 2010, sex, ['prevalence', 'remission'])
    save_country_level_posterior(dm, region, year, sex, ['prevalence', 'remission'])  #'prevalence incidence remission excess-mortality duration mortality relative-risk'.split())


    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])
    dm.save('dm-%d-posterior-%s-%s-%s.json' % (dm.id, region, sex, year), keys_to_save=keys)

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

    dm = fit_without_confrontation(id, options.region, options.sex, options.year)
    return dm

if __name__ == '__main__':
    dm = main()

