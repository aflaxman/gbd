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


    ## separate out prevalence data
    prev_data = [d for d in dm.data if dismod3.gbd_disease_model.relevant_to(d, 'prevalence', region, year, sex)]
    dm.data = [d for d in dm.data if not d in prev_data]


    ### setup the generic disease model (without prevalence data)
    import dismod3.gbd_disease_model as model
    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])
    dm.calc_effective_sample_size(dm.data)
    dm.vars = model.setup(dm, keys)


    ## override the birth prevalence prior, based on the withheld prevalence data
    logit_C_0 = dm.vars[dismod3.utils.gbd_key_for('bins', region, year, sex)]['initial'][2]
    assert len(prev_data) == 1, 'should be a single prevalance datum'
    d = prev_data[0]

    mu_logit_C_0 = mc.logit(dm.value_per_1(d))
    lb, ub = dm.bounds_per_1(d)
    sigma_logit_C_0 = (mc.logit(ub) - mc.logit(lb)) / (2 * 1.96)
    print 'mu_logit_C_0_post:', mu_logit_C_0
    print 'ui_logit_C_0_post:', mc.logit(lb), mc.logit(ub)

    logit_C_0.parents['mu'] = mu_logit_C_0
    logit_C_0.parents['tau'] = sigma_logit_C_0**-2


    ### fit the model
    ## first generate decent initial conditions
    model.fit(dm, method='map', keys=keys, verbose=1)

    ## then sample the posterior via MCMC
    model.fit(dm, method='mcmc', keys=keys, iter=10000, thin=5, burn=5000, verbose=1, dbname='/dev/null')

    print 'mu_logit_C_0_post:', logit_C_0.stats()['mean']
    print 'ui_logit_C_0_post:', logit_C_0.stats()['95% HPD interval']

    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])
    dm.save('dm-%d-posterior-%s-%s-%s.json' % (id, region, sex, year), keys_to_save=keys)

    return dm

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

