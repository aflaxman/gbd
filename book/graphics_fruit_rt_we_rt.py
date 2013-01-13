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

# set font
book_graphics.set_font()

def my_axis(ymax):
    pl.axis([-5,105,-ymax/10.,ymax])
   
def load_USA_new_model():
    # example of predicting out-of-sample with a ln_ASDR covariate
    try:
        model = dismod3.data.load('/home/j/Project/dismod/output/dm-35020')
    except:
        model = dismod3.data.load('/home/j/Project/dismod/dismod_status/prod/dm-35020')
    model.keep(areas=['USA'], sexes=['male'], start_year=1997)
    model.input_data = model.input_data.drop(['z_cv_natl_rep','x_cv_diet_assess_method','x_cv_met_suboptimal','x_cv_natl_rep','x_fao_factor1','x_fao_factor2','x_fao_factor4','x_ln_LDI_pc','x_ln_fruits'], 1)
    return model

def load_we_new_model():
    # example of predicting out-of-sample with a ln_ASDR covariate
    try:
        model = dismod3.data.load('/home/j/Project/dismod/output/dm-35020')
    except:
        model = dismod3.data.load('/home/j/Project/dismod/dismod_status/prod/dm-35020')
    model.keep(areas=['europe_western'], sexes=['male'], start_year=1997)
    model.input_data = model.input_data.drop(['z_cv_natl_rep','x_cv_diet_assess_method','x_cv_met_suboptimal','x_cv_natl_rep','x_fao_factor1','x_fao_factor2','x_fao_factor4','x_ln_LDI_pc','x_ln_fruits'], 1)
    return model

best_model = load_USA_new_model()
we_model = load_we_new_model()

output = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_USA_we.csv')

# figure 

pl.figure(**book_graphics.full_page_params)
pl.subplot(1,2,1)
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
pl.legend(loc='upper center', bbox_to_anchor=(.2,-.23), fancybox=True, shadow=True) 

pl.subplot(1,2,2)
pl.plot(pl.array(output['USA_model'])-pl.array(output['USA_log_model']), 'k--', linewidth=3, label='Negative-binomial -\n lognormal posterior mean')
pl.plot(pl.array(output['USA_model'])-pl.array(output['USA_norm_model']), 'k:', linewidth=3, label='Negative-binomial -\n normal posterior mean')

pl.yticks([-.0001, -.0, .0001, .0002])
pl.axis([-5,105,-.00015,.00025])

pl.subplots_adjust(wspace=.5, top=.99, bottom=.62)
pl.legend(loc='upper center', bbox_to_anchor=(.2,-.23), fancybox=True, shadow=True) 

pl.savefig('book/graphics/fruit-rate_type.pdf')
pl.savefig('book/graphics/fruit-rate_type.png')

