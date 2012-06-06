""" Validate Consistent Model based on Major Depression dataset
"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import networkx as nx
import pymc as mc
import pandas

import data
import consistent_model
import fit_model
import graphics
import data_simulation
import covariate_model
import dismod3

reload(consistent_model)
reload(fit_model)
reload(data_simulation)
reload(graphics)

pl.seterr('ignore')

types_to_resample = 'i r f p pf m_with rr'.split()


def resample(data):
    if len(data) == 0:
        return data

    delta_true = .1
    p = data['mu_pred']+1.e-6


    # TODO: abstract this block of code into rate_model.py; it is also called in data_model.py
    ## ensure that all data has uncertainty quantified appropriately
    # first replace all missing se from ci
    missing_se = pl.isnan(data['standard_error']) | (data['standard_error'] <= 0)
    data['standard_error'][missing_se] = (data['upper_ci'][missing_se] - data['lower_ci'][missing_se]) / (2*1.96)

    # then replace all missing ess with se
    missing_ess = pl.isnan(data['effective_sample_size'])
    data['effective_sample_size'][missing_ess] = data['value'][missing_ess]*(1-data['value'][missing_ess])/data['standard_error'][missing_ess]**2

    # warn and drop data that doesn't have effective sample size quantified, or is is non-positive
    missing_ess = pl.isnan(data['effective_sample_size']) | (data['effective_sample_size'] < 0)
    if sum(missing_ess) > 0:
        print 'WARNING: %d rows of data has invalid quantification of uncertainty.' % sum(missing_ess)
        data['effective_sample_size'][missing_ess] = 1.0

    n = data['effective_sample_size']

    data['true'] = p
    data['value'] = (1.0 * mc.rnegative_binomial(n*p, delta_true*n*p)) / n

    # uncomment below to test the effect of having very wrong data
    #data['value'] = 0.
    #data['effective_sample_size'] = 1.e6

    return data
    


def simulate_data(dm, area, sex, year):
    """ generate simulated data, based on mu_pred of dm.vars[t]['data']
    """
    param_type = dict(i='incidence', p='prevalence', r='remission', f='excess-mortality', rr='relative-risk', pf='prevalence_x_excess-mortality', m_with='mortality')


    # TODO: abstract this block of code into data.py, it is also called in fit_posterior.py
    # select data that is about areas in this region, recent years, and sex of male or total only
    predict_area = area
    predict_year = year
    predict_sex = sex
    model = dm.model
    assert predict_area in model.hierarchy, 'region %s not found in area hierarchy' % predict_area
    subtree = nx.traversal.bfs_tree(model.hierarchy, predict_area)
    relevant_rows = [i for i, r in model.input_data.T.iteritems() \
                         if (r['area'] in subtree or r['area'] == 'all')\
                         and ((predict_year >= 1997 and r['year_end'] >= 1997) or
                              (predict_year <= 1997 and r['year_start'] <= 1997)) \
                         and r['sex'] in [predict_sex, 'total']]
    model.input_data = model.input_data.ix[relevant_rows]

    # replace area 'all' with predict_area
    model.input_data['area'][model.input_data['area'] == 'all'] = predict_area


    input_data = pandas.DataFrame()
    for t in types_to_resample:
        data_t = model.get_data(t)
        if t == 'rr':
            input_data = input_data.append(data_t)  # TODO: resample rr, smr, X with appropriate noise model
        else:
            input_data = input_data.append(resample(data_t))
    # copy m_all data over
    input_data = input_data.append(dm.model.input_data[dm.model.input_data['data_type']=='m_all'])
    dm.model.input_data = input_data


    true = {}
    for t in types_to_resample:
        key = dismod3.utils.gbd_key_for(param_type[t], area, year, sex)
        true[t] = dm.get_mcmc('mean', key)
    dm.true = true


    emp_priors = {}
    for t in 'i r f p rr pf'.split():
        key = dismod3.utils.gbd_key_for(param_type[t], area, year, sex)
        mu = dm.get_mcmc('emp_prior_mean', key)
        sigma = dm.get_mcmc('emp_prior_std', key)
        
        if len(mu) == 101 and len(sigma) == 101:
            emp_priors[t, 'mu'] = mu
            emp_priors[t, 'sigma'] = sigma
    dm.emp_priors = emp_priors
    dm.vars = consistent_model.consistent_model(dm.model, area, sex, year, emp_priors)

    return dm


def fit_simulated(dm, area, sex, year):
    #dm.map, dm.mcmc = fit_model.fit_consistent_model(dm.vars, iter=101, burn=0, thin=1, tune_interval=100)
    dm.map, dm.mcmc = fit_model.fit_consistent_model(dm.vars, iter=10000, burn=5000, thin=25, tune_interval=100)

    posteriors = {}
    for t in 'i r f p rr pf'.split():
        est_k = covariate_model.predict_for(dm.model, area, sex, year, area, sex, year, 1., dm.vars[t], 0., pl.inf)
        posteriors[t] = est_k
    dm.posteriors = posteriors

def store_results(dm, area, sex, year):
    types_to_plot = 'p i r rr'.split()

    graphics.plot_convergence_diag(dm.vars)
    pl.clf()
    for i, t in enumerate(types_to_plot):
        pl.subplot(len(types_to_plot), 1, i+1)
        graphics.plot_data_bars(dm.model.get_data(t))
        pl.plot(range(101), dm.emp_priors[t, 'mu'], linestyle='dashed', color='grey', label='Emp. Prior', linewidth=3)
        pl.plot(range(101), dm.true[t], 'b-', label='Truth', linewidth=3)
        pl.plot(range(101), dm.posteriors[t].mean(0), 'r-', label='Estimate', linewidth=3)

        pl.errorbar(range(101), dm.posteriors[t].mean(0), yerr=1.96*dm.posteriors[t].std(0), fmt='r-', linewidth=1, capsize=0)

        pl.ylabel(t)
        graphics.expand_axis()
    
    pl.legend(loc=(0.,-.95), fancybox=True, shadow=True)
    pl.subplots_adjust(hspace=0, left=.1, right=.95, bottom=.2, top=.95)
    pl.xlabel('Age (Years)')
    pl.show()

    model = dm
    model.mu = pandas.DataFrame()
    for t in types_to_plot:
        model.mu = model.mu.append(pandas.DataFrame(dict(true=dm.true[t],
                                                         mu_pred=dm.posteriors[t].mean(0),
                                                         sigma_pred=dm.posteriors[t].std(0))),
                                   ignore_index=True)
    data_simulation.add_quality_metrics(model.mu)
    print '\nparam prediction bias: %.5f, MARE: %.3f, coverage: %.2f' % (model.mu['abs_err'].mean(),
                                                                         pl.median(pl.absolute(model.mu['rel_err'].dropna())),
                                                                         model.mu['covered?'].mean())
    print

    data_simulation.initialize_results(model)
    data_simulation.add_to_results(model, 'mu')
    data_simulation.finalize_results(model)

    print model.results

    return model


if __name__ == '__main__':
    region, sex, year = 'north_america_high_income', 'male', 1990

    import fit_posterior, upload_fits
    import data
    import simplejson as json

    ## load the model from disk or from web
    dm = dismod3.load_disease_model(24842)
    dm.model = data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))
    data = upload_fits.merge_data_csvs(24842)
    dm.model.input_data['mu_pred'] = data['mu_pred']

    simulate_data(dm, region, sex, year)
    fit_simulated(dm, region, sex, year)
    store_results(dm, region, sex, year)
