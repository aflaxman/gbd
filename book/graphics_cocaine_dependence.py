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
    pl.axis([11,69,-ymax/10.,ymax])
    pl.grid()
	
def subtitle(s):
    """ title where the panel names appear within each panel"""
    l,r,b,t=pl.axis()
    x = l + (r-l)*.05
    y = t - (t-b)*.05
    pl.text(x, y, s, ha='left', va='top')

def load_new_model():
    orig_model = dismod3.data.load('/home/j/Project/dismod/output/dm-40358') 
    
    #orig_model.input_data['area']='all'
    orig_model.keep(areas=['USA'])
    
    #orig_model.input_data = orig_model.input_data[orig_model.input_data['x_cv_canuse'] == 0]
    #orig_model.input_data = orig_model.input_data[orig_model.input_data['x_cv_past_year'] == 1]
    #orig_model.input_data = orig_model.input_data.drop(['x_cv_canuse', 'x_cv_past_year'], axis=1)
    
    # remove smoothness parameter initially
    orig_model.parameters['p']['smoothness']['amount'] = 'Slightly'
    orig_model.parameters['p']['random_effects'] = {}
    orig_model.parameters['p']['parameter_age_mesh'] = [15,20,25,30,40,50,65]

    return orig_model

best_model = load_new_model()

output = pandas.read_csv('/homes/peterhm/gbd/book/applications-data_cocaine_dep.csv')

# figure cocaine_dependence-data
pl.figure(**book_graphics.full_page_params)

dismod3.graphics.plot_data_bars(best_model.get_data('p'))
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .005, .01, .015], [0, .5, 1.0, 1.5])
my_axis(.016)
pl.savefig('/homes/peterhm/gbd/book/applications/cocaine_dependence-data.pdf')
pl.savefig('/homes/peterhm/gbd/book/applications/cocaine_dependence-data.png')

# figure cocaine_dependence-knots
pl.figure(**book_graphics.full_page_params)

x = best_model.parameters['p']['parameter_age_mesh']

dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')

pl.plot(x, pl.array(output['best_model'])[x], 'k-', label='7 knots', linewidth=3)
pl.plot(x, pl.array(output['best_model_l'])[x], 'k-', label='7 knots 95% \nHPD interval', linewidth=1)
pl.plot(x, pl.array(output['best_model_u'])[x], 'k-', linewidth=1)
pl.plot(x, pl.array(output['less_model'])[x], 'k:', label='4 knots', linewidth=3)
pl.plot(x, pl.array(output['more_model'])[x], 'k--', label='11 knots', linewidth=3)
#plot(most_model.vars['p']['mu_age'].stats()['quantiles'][50], 'k-.', label='Most knots')
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .0025, .005, .0075, .010, 0.0125, .015], [0, .25, .5, .75, 1.0, 1.25, 1.5])
my_axis(.017)
pl.legend(loc='upper right', fancybox=True, shadow=True)

pl.savefig('/homes/peterhm/gbd/book/applications/cocaine_dependence-knots.pdf')
pl.savefig('/homes/peterhm/gbd/book/applications/cocaine_dependence-knots.png')

# figure cocaine_dependence-smoothing
more_model_slightly = load_new_model()
more_model_slightly.parameters['p']['parameter_age_mesh'] = range(15,66,5)

pl.figure(**book_graphics.full_page_params)
x = more_model_slightly.parameters['p']['parameter_age_mesh']

dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')

pl.plot(x, pl.array(output['none'])[x], 'k:', label='$\sigma$ = $\infty$', linewidth=3)
pl.plot(x, pl.array(output['moderately'])[x], 'k-.', label='$\sigma$ = 0.05', linewidth=3)
pl.plot(x, pl.array(output['very'])[x], 'k--', label='$\sigma$ = 0.005', linewidth=3)

pl.plot(x, pl.array(output['slightly'])[x], 'k-', label='$\sigma$ = 0.1', linewidth=3)
pl.plot(x, pl.array(output['slightly_l'])[x], 'k-', label='$\sigma$ = 0.1 95% \nHPD interval', linewidth=1)
pl.plot(x, pl.array(output['slightly_u'])[x], 'k-', linewidth=1)

pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .004, .008, .012, .016], [0, 0.4, 0.8, 1.2, 1.6])
my_axis(.018)
pl.legend(loc='upper right', fancybox=True, shadow=True, ncol=2)

pl.savefig('/homes/peterhm/gbd/book/applications/cocaine_dependence-smoothing.pdf')
pl.savefig('/homes/peterhm/gbd/book/applications/cocaine_dependence-smoothing.png')

