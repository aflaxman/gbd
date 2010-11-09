""" Script to fit the North America High Income Hep C data for 1990

"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 


from pylab import *
import pymc as mc
    
from dismod3.disease_json import DiseaseJson
from dismod3 import neg_binom_model
import dismod3.utils
import test_model

def hep_c_fit():
    """ Test fit for prevalence only"""

    # load model to fit
    dm = DiseaseJson(file('tests/hep_c.json').read())

    # adjust the expert priors
    dm.params['global_priors']['heterogeneity']['prevalence'] = 'Very'
    # TODO: model adjusting other covariates
    
    # TODO: consider covariates

    # select relevant prevalence data
    # TODO: streamline data selection functions
    dm.data = [d for d in dm.data if
               dismod3.utils.clean(d['gbd_region']) == 'north_america_high_income'
               and float(d['year_start']) <= 1997]

    # create, fit, and save rate models
    dm.vars = {}

    keys = dismod3.utils.gbd_keys(type_list=['prevalence'],
                                  region_list=['north_america_high_income'],
                                  year_list=[1990],)
    for k in keys:
        # TODO: consider how to do this for models that use the complete disease model
        dm.vars[k] = neg_binom_model.setup(dm, k, dm.data)

        # TODO: consider adding hierarchical similarity priors for the male and female models
        dm.mcmc = mc.MCMC(dm.vars)
        dm.mcmc.sample(iter=20000, burn=10000, thin=10, verbose=1)

        # make map object so we can compute AIC and BIC
        dm.map = mc.MAP(dm.vars)
        dm.map.fit()

        # save the results in the disease model
        neg_binom_model.store_mcmc_fit(dm, k, dm.vars[k])

        # check autocorrelation to confirm chain has mixed
        test_model.summarize_acorr(dm.vars[k]['rate_stoch'].trace())

        # generate plots of results
        dismod3.tile_plot_disease_model(dm, [k], defaults={'ymax':.15})
        dm.savefig('dm-%d-posterior-%s.%f.png' % (dm.id, k, random()))

        # summarize fit quality graphically, as well as parameter posteriors
        dismod3.plotting.plot_posterior_predicted_checks(dm, k)
        dm.savefig('dm-%d-check-%s.%f.png' % (dm.id, k, random()))
    return dm

if __name__ == '__main__':
    dm = hep_c_fit()
    
