import sys
sys.path += ['../gbd', '../gbd/book', '../dm3-computation_only/', '../dm3-computation_only/book']
import pylab as pl
import pymc as mc
import pandas

import dismod3
reload(dismod3)

import book_graphics
reload(book_graphics)
import matplotlib as mpl

# make all fonts bigger, etc

mpl.rcParams['axes.titlesize'] = 'xx-large'
mpl.rcParams['axes.labelsize'] = 'xx-large'

mpl.rcParams['xtick.labelsize'] = 'x-large'
mpl.rcParams['ytick.labelsize'] = 'x-large'

mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.fontsize'] = 'large'

mpl.rcParams['text.fontsize'] = 12

def my_axis(ymax):
    pl.axis([-5,105,-ymax/10.,ymax])
	
def subtitle(s):
    """ title where the panel names appear within each panel"""
    l,r,b,t=pl.axis()
    x = l + (r-l)*.05
    y = t - (t-b)*.05
    pl.text(x, y, s, ha='left', va='top')

def my_plot_data_bars(df, color, label, style='book'):
    """ Plot some data bars
    Input
    -----
    df : pandas.DataFrame with columns age_start, age_end, value
    """
    data_bars = zip(df['age_start'], df['age_end'], df['value'])

    # show at most 500 bars, to keep things fast
    # TODO: make 500 into an option
    if len(data_bars) > 500:
        import random
        data_bars = random.sample(data_bars, 500)

    # make lists of x and y points, faster than ploting each bar
    # individually
    x = []
    y = []
    for a_0i, a_1i, p_i in data_bars:
        x += [a_0i, a_1i, pl.nan]
        y += [p_i, p_i, pl.nan]

    pl.plot(x, y, 's-', mew=1, mec='w', ms=4, color=color, label=label)
    
def load_USA_new_model():
    # example of predicting out-of-sample with a ln_ASDR covariate
    model = dismod3.data.load('/home/j/Project/dismod/output/dm-35020')
    model.keep(areas=['USA'], sexes=['male'], start_year=1997)
    model.input_data = model.input_data.drop(['z_cv_natl_rep','x_cv_diet_assess_method','x_cv_met_suboptimal','x_cv_natl_rep','x_fao_factor1','x_fao_factor2','x_fao_factor4','x_ln_LDI_pc','x_ln_fruits'], 1)
    return model

def load_we_new_model():
    # example of predicting out-of-sample with a ln_ASDR covariate
    model = dismod3.data.load('/home/j/Project/dismod/output/dm-35020')
    model.keep(areas=['europe_western'], sexes=['male'], start_year=1997)
    model.input_data = model.input_data.drop(['z_cv_natl_rep','x_cv_diet_assess_method','x_cv_met_suboptimal','x_cv_natl_rep','x_fao_factor1','x_fao_factor2','x_fao_factor4','x_ln_LDI_pc','x_ln_fruits'], 1)
    return model

best_model = load_USA_new_model()
we_model = load_we_new_model()

output = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_USA_we.csv')

# figure 

pl.figure(**book_graphics.full_page_params)
dismod3.graphics.plot_data_bars(best_model.get_data('r'), color='grey') 
   
knots = best_model.parameters['r']['parameter_age_mesh']

pl.plot(pl.array(output['USA_model']), 'k-', linewidth=3, label='Negative-binomial posterior mean')
pl.plot(knots, pl.array(output['USA_model_l'])[knots], 'k-', linewidth=1, label='Negative-binomial 95% HPD interval')
pl.plot(knots, pl.array(output['USA_model_u'])[knots], 'k-', linewidth=1)

pl.plot(pl.array(output['USA_log_model']), 'k--', linewidth=3, label='Lognormal posterior mean')
pl.plot(knots, pl.array(output['USA_log_model_l'])[knots], 'k--', linewidth=1, label='Lognormal 95% HPD interval')
pl.plot(knots, pl.array(output['USA_log_model_u'])[knots], 'k--', linewidth=1)

pl.plot(pl.array(output['USA_norm_model']), 'k:', linewidth=3, label='Normal posterior mean')
pl.plot(knots, pl.array(output['USA_norm_model_l'])[knots], 'k:', linewidth=1, label='Normal 95% HPD interval')
pl.plot(knots, pl.array(output['USA_norm_model_u'])[knots], 'k:', linewidth=1)

pl.xlabel('Age (years)')
pl.ylabel('Consumption (kg/d)')
pl.yticks([0, .006, .012, .018, .024], [0, 0.06, 0.12, 0.18, 0.24])
my_axis(.026)
pl.legend(loc='lower right', fancybox=True, shadow=True)
pl.grid()

pl.savefig('/homes/peterhm/gbd/book/applications/fruit-rate_type.pdf')
pl.savefig('/homes/peterhm/gbd/book/applications/fruit-rate_type.png')

# figure-we_rate_type
isl_gbr = load_we_new_model()
isl_gbr.keep(areas=['ISL', 'GBR'])

pl.figure(**book_graphics.full_page_params)
my_plot_data_bars(we_model.get_data('r'), color='grey', label='Western Europe') 
my_plot_data_bars(isl_gbr.get_data('r'), color='black', label='Great Britain and Iceland')

knots = we_model.parameters['r']['parameter_age_mesh']

pl.plot(pl.array(output['we_model']), 'k-', linewidth=3, label='Negative-binomial posterior mean')
pl.plot(knots, pl.array(output['we_model_l'])[knots], 'k-', linewidth=1, label='Negative-binomial 95% HPD interval')
pl.plot(knots, pl.array(output['we_model_u'])[knots], 'k-', linewidth=1)

pl.plot(pl.array(output['we_log_model']), 'k--', linewidth=3, label='Lognormal posterior mean')
pl.plot(knots, pl.array(output['we_log_model_l'])[knots], 'k--', linewidth=1, label='Lognormal 95% HPD interval')
pl.plot(knots, pl.array(output['we_log_model_u'])[knots], 'k--', linewidth=1)

pl.plot(pl.array(output['we_norm_model']), 'k:', linewidth=3, label='Normal posterior mean')
pl.plot(knots, pl.array(output['we_norm_model_l'])[knots], 'k:', linewidth=1, label='Normal 95% HPD interval')
pl.plot(knots, pl.array(output['we_norm_model_u'])[knots], 'k:', linewidth=1)

pl.xlabel('Age (years)')
pl.ylabel('Consumption (kg/d)')
pl.yticks([0, .01, .02, .03, .04], [0, 0.1, 0.2, 0.3, 0.4])
my_axis(.06)
pl.legend(loc='upper right', fancybox=True, shadow=True)
pl.grid()

pl.subplots_adjust(hspace=.3)

pl.savefig('/homes/peterhm/gbd/book/applications/fruit-we_rate_type.pdf')
pl.savefig('/homes/peterhm/gbd/book/applications/fruit-we_rate_type.png')