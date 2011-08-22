#!/usr/bin/python2.5
""" Generate empirical prior of specified parameter type

Expects the disase model json to be saved already.
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

import dismod3
import Matplot

import pymc as mc
import pylab as pl

def fit_world(dm):
    """ Fit consistent for all data in world

    Parameters
    ----------
    id : int
      The model id number for the job to fit

    Example
    -------
    >>> import fit_world
    >>> dm = fit_world.dismod3.load_disease_model(1234)
    >>> fit_world.fit_world(dm)
    """
    region = 'world'
    year = 'all'
    sex = 'all'

    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])

    print 'initializing model vars... ',
    dm.calc_effective_sample_size(dm.data)
    for k in keys:
        dm.fit_initial_estimate(k)

    dm.vars = dismod3.generic_disease_model.setup(dm, dismod3.utils.gbd_key_for(type='%s', region=region, year=year, sex=sex))

    verbose=1
    def map_fit(stoch_names):
        print '\nfitting', ' '.join(stoch_names)
        key = dismod3.utils.gbd_key_for('%s', region, year, sex)

        # fit all parameters besides over-dispersion
        import pymc as mc
        map = mc.MAP([dm.vars[key%type][subkey] for type in stoch_names for subkey in dm.vars[key%type] if not subkey.startswith('log_dispersion')])
        try:
                                                      map.fit(method='fmin_powell', verbose=verbose)
        except KeyboardInterrupt:
            print 'User halted optimization routine before optimal value found'

        for type in stoch_names:
            for subkey in ['age_coeffs_mesh', 'dispersion']:
                if subkey in dm.vars[key%type]:
                    print key%type, subkey, pl.atleast_1d(dm.vars[key%type][subkey].value).round(0)

        return map

    # use map as initial values
    print 'initializing MAP object... ',
    dm.map = map_fit('incidence remission excess-mortality mortality relative-risk smr prevalence_x_excess-mortality duration bins prevalence'.split())
    print 'initialization completed'

    # fit with mcmc
    import pymc as mc
    dm.mcmc = mc.MCMC(dm.vars)
    for k in keys:
        if 'dispersion_step_sd' in dm.vars[k]:
            dm.mcmc.use_step_method(mc.Metropolis, dm.vars[k]['log_dispersion'],
                                    proposal_sd=dm.vars[k]['dispersion_step_sd'])
        if 'age_coeffs_mesh_step_cov' in dm.vars[k]:
            dm.mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['age_coeffs_mesh'],
                                    cov=dm.vars[k]['age_coeffs_mesh_step_cov'], verbose=0)

            # TODO: make a wrapper function for handling this adaptive metropolis setup
            stoch_list = [dm.vars[k]['study_coeffs'], dm.vars[k]['region_coeffs'], dm.vars[k]['age_coeffs_mesh']]
            d1 = len(dm.vars[k]['study_coeffs'].value)
            d2 = len(dm.vars[k]['region_coeffs_step_cov'])
            d3 = len(dm.vars[k]['age_coeffs_mesh_step_cov'])
            C = pl.eye(d1+d2+d3)
            C[d1:(d1+d2), d1:(d1+d2)] = dm.vars[k]['region_coeffs_step_cov']
            C[(d1+d2):(d1+d2+d3), (d1+d2):(d1+d2+d3)] = dm.vars[k]['age_coeffs_mesh_step_cov']
            C *= .01
            dm.mcmc.use_step_method(mc.AdaptiveMetropolis, stoch_list, cov=C)

            # more step methods
            dm.mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['age_coeffs_mesh'], cov=dm.vars[k]['age_coeffs_mesh_step_cov'])
        if isinstance(dm.vars[k].get('study_coeffs'), mc.Stochastic):
            dm.mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['study_coeffs'])
        if isinstance(dm.vars[k].get('region_coeffs'), mc.Stochastic):
            dm.mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['region_coeffs'], cov=dm.vars[k]['region_coeffs_step_cov'])
    dm.mcmc.sample(iter=20000, burn=10000, thin=10, verbose=verbose)

    # generate plots
    for k in dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex]):
        if dm.vars[k].get('data'):
            dismod3.plotting.plot_posterior_predicted_checks(dm, k)
            dm.savefig('dm-%d-check-%s.png' % (dm.id, k))

    dir = dismod3.settings.JOB_WORKING_DIR % dm.id
    for k in keys:
        t,r,y,s = dismod3.utils.type_region_year_sex_from_key(k)
        if t in ['incidence', 'prevalence', 'remission', 'excess-mortality']:
            for s in 'dispersion age_coeffs_mesh study_coeffs region_coeffs'.split():
                if s in dm.vars[k] and isinstance(dm.vars[k][s], mc.Node):
                    try:
                        Matplot.plot(dm.vars[k][s], path='%s/image/mcmc_diagnostics/'%dir, common_scale=False)
                        pass
                    except Exception, e:
                        print e

    # save the results
    for param_type in 'prevalence incidence remission excess-mortality'.split():
        key = dismod3.utils.gbd_key_for(param_type, region, year, sex)
        prior_vals = dict(
            alpha=list(dm.vars[key]['region_coeffs'].stats()['mean']),
            beta=list(pl.atleast_1d(dm.vars[key]['study_coeffs'].stats()['mean'])),
            gamma=list(dm.vars[key]['age_coeffs'].stats()['mean']),
            delta=float(dm.vars[key]['dispersion'].stats()['mean']))

        prior_vals.update(
            sigma_alpha=list(dm.vars[key]['region_coeffs'].stats()['standard deviation']),
            sigma_beta=list(pl.atleast_1d(dm.vars[key]['study_coeffs'].stats()['standard deviation'])),
            sigma_gamma=list(dm.vars[key]['age_coeffs'].stats()['standard deviation']),
            sigma_delta=float(dm.vars[key]['dispersion'].stats()['standard deviation']))
        dm.set_empirical_prior(param_type, prior_vals)

    for param_type in 'prevalence incidence remission excess-mortality'.split():
        key = dismod3.utils.gbd_key_for(param_type, region, year, sex)
        trace = zip(dm.vars[key]['region_coeffs'].trace(),
                    dm.vars[key]['study_coeffs'].trace(),
                    dm.vars[key]['age_coeffs'].trace())
        bound_func = dm.vars[key]['bounds_func']
        for r in dismod3.gbd_regions:
            print 'predicting rates for %s' % r
            for y in dismod3.gbd_years:
                for s in dismod3.gbd_sexes:
                    pred_key = dismod3.utils.gbd_key_for(param_type, r, y, s)
                    rate_trace = []
                    for a, b, g in trace:
                        rate_trace.append(dismod3.neg_binom_model.predict_region_rate(pred_key,
                                                                                      alpha=a,
                                                                                      beta=b,
                                                                                      gamma=g,
                                                                                      covariates_dict=dm.get_covariates(),
                                                                                      derived_covariate=dm.get_derived_covariate_values(),
                                                                                      bounds_func=bound_func,
                                                                                      ages=dm.get_estimate_age_mesh()))
                    mu = pl.mean(rate_trace, axis=0)
                    dm.set_mcmc('emp_prior_mean', pred_key, mu)
    
    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    dm.save('dm-%d-prior.json' % dm.id)

def main():
    import optparse

    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    dm = dismod3.load_disease_model(id)
    fit_world(dm)
    return dm
      

if __name__ == '__main__':
    dm = main()
