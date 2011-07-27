#!/usr/bin/python2.5
""" Generate a posterior estimate for a specific region, sex, and year
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

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

    dm = dismod3.load_disease_model(id)
    
    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])

    # fit the model
    dir = dismod3.settings.JOB_WORKING_DIR % id
    import dismod3.gbd_disease_model as model
    model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    ## then sample the posterior via MCMC
    model.fit(dm, method='mcmc', keys=keys, iter=100000, thin=50, burn=50000, verbose=1, dbname='/dev/null')

    # generate plots of results
    dismod3.tile_plot_disease_model(dm, keys, defaults={})
    # TODO: refactor naming into its own function (disease_json.save_image perhaps)
    dm.savefig('dm-%d-posterior-%s.png' % (id, dismod3.utils.gbd_key_for('all', region, year, sex)))

    # summarize fit quality graphically, as well as parameter posteriors
    for k in dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex]):
        if dm.vars[k].get('data'):
            dismod3.plotting.plot_posterior_predicted_checks(dm, k)
            dm.savefig('dm-%d-check-%s.png' % (dm.id, k))

    # make a rate_type_list
    rate_type_list = ['incidence', 'prevalence', 'remission', 'excess-mortality', 'mortality', 'duration']
    rate_type_list = ['prevalence']

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

