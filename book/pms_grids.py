import sys
sys.path += ['..']

faster_run_flag = False

import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}

### @export 'data'
grids = dict(a=[0,15,20,40,50,100],
             b=[0,15,17,40,50,100],
             c=[0] + range(15,51,5) + [100],
             d=[0] + range(15,51,2) + [100])

region = 'north_america_high_income'
sex = 'female'
year = '2005'

for grid in grids.keys():
    # load model
    dm = dismod3.load_disease_model(16370)

    # set expert priors and other model parameters
    dm.set_param_age_mesh(grids[grid])

    dm.params['global_priors']['level_value']['incidence']['age_before'] = 15
    dm.params['global_priors']['level_value']['incidence']['age_after'] = 50
    dm.params['global_priors']['smoothness']['incidence']['age_start'] = 15

    dm.params['global_priors']['level_value']['remission']['age_before'] = 40
    dm.params['global_priors']['level_bounds']['remission']['upper'] = 10.
    
    dm.params['global_priors']['level_value']['excess_mortality']['age_before'] = 101
    
    dm.params['global_priors']['level_value']['prevalence']['age_before'] = 15
    dm.params['global_priors']['level_value']['prevalence']['age_after'] = 50
    dm.params['global_priors']['smoothness']['prevalence']['age_start'] = 15
    
    dm.params['covariates']['Country_level']['LDI_id_Updated_7July2011']['rate']['value'] = 0


    # clear any fit and priors
    dm.clear_fit()
    dm.clear_empirical_prior()
    dismod3.neg_binom_model.covariate_hash = {}

    import fit_posterior
    data = dm.data
    fit_posterior.fit_posterior(dm, region, sex, year, map_only=faster_run_flag, store_results=False)
    dm.data = data # put data back in
    results[grid] = dm
    try:
        results['dic_%s'%grid] = dm.mcmc.dic
    except Exception, e:
        print e
        results['dic_%s'%grid] = 'TK'


pl.figure(**book_graphics.quarter_page_params)
for ii, grid in enumerate('abcd'):
    dm = results[grid]
    pl.subplot(1,4,ii+1)
    dismod3.plotting.plot_intervals(dm, [d for d in dm.data if dm.relevant_to(d, 'prevalence', region, year, sex)],
                                    color='black', print_sample_size=False, alpha=1., plot_error_bars=False,
                                    linewidth=2)
    book_graphics.plot_rate(dm, dismod3.utils.gbd_key_for('prevalence', region, year, sex), linestyle='-')
    pl.axis([10, 60, -.05, .8])
    pl.xlabel('Age (Years)')
    pl.xticks([15,35,55])
    if ii == 0:
        pl.yticks([0, .2, .4, .6], [0, 20, 40, 60])
        pl.ylabel('Prevalence (per 100)')
    else:
        pl.yticks([0, .2, .4, .6], ['', '', '', '']) 
    pl.text(12, .75, '(%s)'% grid, va='top', ha='left')


pl.subplots_adjust(wspace=0, bottom=.15, left=.1, right=.99, top=.97)
pl.savefig('pms_grids.pdf')
book_graphics.save_json('pms_grids.json', vars())
