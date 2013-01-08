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

def load_new_model():
    # example of predicting out-of-sample with a ln_ASDR covariate
    model = dismod3.data.load('/home/j/Project/dismod/output/dm-37005')
    #model.input_data = model.input_data.drop(['x_health_system_access'], axis=1)
    return model

# figure cirrhosis-data
best_model = load_new_model()

pl.figure(**book_graphics.half_page_params)

pl.subplot(1,2,1)
dismod3.graphics.plot_data_bars(best_model.get_data('p'))
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .001, .002, .003, .004], [0, 0.1, 0.2, 0.3, 0.4])
my_axis(.0045)
subtitle('(a)')
pl.grid()

pl.subplot(1,2,2)
dismod3.graphics.plot_data_bars(best_model.get_data('pf'))
pl.xlabel('Age (years)')
pl.ylabel('Cause-specific mortality \n (per 10,000 PY)'+'\n\n', ha='center')
pl.yticks([0, .0008, .0016, .0024, .0032], [0, 8, 16, 24, 32])
my_axis(.0035)
pl.subplots_adjust(hspace=.35)
pl.subplots_adjust(wspace=.35)
subtitle('(b)')
pl.grid()

pl.subplots_adjust(wspace=.35, hspace=.35, bottom=.14)

pl.savefig('/homes/peterhm/gbd/book/applications/cirrhosis-data.pdf')
pl.savefig('/homes/peterhm/gbd/book/applications/cirrhosis-data.png')

# figure cirrhosis-lnASDR_v_prev
pl.figure(**book_graphics.full_page_params)

d=best_model.get_data('p')
d = d[d['age_start'] == 60]
x = d['x_lnASDR_B05'].__array__()
y = d['value'].__array__()
pl.plot(pl.exp(x), y, 'ks')
pl.xlabel('Age-standardized death rate (per 100,000 PY)')
pl.ylabel('Prevalence (%)')
pl.xticks([0, .00002, .00004, .00006, .00008], [0, 2, 4, 6, 8])
pl.yticks([0, .001, .002, .003, .004], [0, 0.1, 0.2, 0.3, 0.4])
pl.axis([-.000001, .00009, -.00005, .0041])
pl.grid()

pl.savefig('/homes/peterhm/gbd/book/applications/cirrhosis-lnASDR_v_prev.pdf')
pl.savefig('/homes/peterhm/gbd/book/applications/cirrhosis-lnASDR_v_prev.png')

# figure cirrhosis-prev_est
output = pandas.read_csv('/homes/peterhm/gbd/book/applications-cirrhosis.csv')

pl.figure(**book_graphics.half_page_params)

param_list = [dict(type='p', title='(a)', ylabel='Prevalence (%)', yticks=([0, .25, .5, .75, 1.], [0, 25, 50, 75, 100]), axis=[-5,105,-.105,1.05], loc='upper left'),
          dict(type='i', title='(b)', ylabel='Incidence (Per 100 PY)', yticks=([0, .1, .2, .3, .4], [0, 10, 20, 30, 40]), axis=[-5,105,-.045,.45], loc='upper left'),
          #dict(type='f', title='(c)', ylabel='Excess mortality (Per 100 PY)', yticks=([0, .002, .004, .006, .008], [0, 2, 4, 6, 8]), axis=[-5,105,-.0018,.009], loc='upper left' ),
          #dict(type='r', title='(c)', ylabel='Excess mortality (Per 100 PY)', yticks=([0, .002, .004, .006, .008], [0, 2, 4, 6, 8]), axis=[-5,105,-.0018,.009], loc='upper left' ),
          #dict(type='m_with', title='(d)', ylabel='With-condition mortality (Per 100 PY)', yticks=([0, .1, .2, .3, .4], [0, 10, 20, 30, 40]), axis=[-5,105,-.045,.45], loc='upper left'),
          ]

x = best_model.parameters['p']['parameter_age_mesh']

pl.subplot(1,2,1)
usa = load_new_model()
usa.keep(areas=['USA'], sexes=['male'])
dismod3.graphics.plot_data_bars(usa.get_data('p'), color='grey') 

pl.plot(pl.arange(101), pl.array(output['USA']), 'k-', linewidth=2, label='Posterior mean')
pl.plot(x, pl.array(output['USA_l'])[x], 'k-', linewidth=1, label='95% HPD interval')
pl.plot(x, pl.array(output['USA_u'])[x], 'k-', linewidth=1)
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .0003, .0006, .0009, .0012], [0, 0.03, 0.06, 0.09, 0.12])
my_axis(.0015)
subtitle('(a)')
pl.grid()

pl.subplot(1,2,2)
egy = load_new_model()
egy.keep(areas=['EGY'], sexes=['male'])
dismod3.graphics.plot_data_bars(egy.get_data('p'), color='grey') 
pl.plot(pl.arange(101), pl.array(output['EGY']), 'k-', linewidth=2, label='Posterior mean')
pl.plot(x, pl.array(output['EGY_l'])[x], 'k-', linewidth=1, label='95% HPD interval')
pl.plot(x, pl.array(output['EGY_u'])[x], 'k-', linewidth=1)
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .002, .004, .006, .008], [0, 0.2, 0.4, 0.6, 0.8])
my_axis(.011)
subtitle('(b)')
pl.grid()
pl.legend(bbox_to_anchor=(.42, 0, .5, .94), bbox_transform=pl.gcf().transFigure, fancybox=True, shadow=True)

pl.subplots_adjust(hspace=.35,wspace=.35,bottom=.14)

pl.savefig('/homes/peterhm/gbd/book/applications/cirrhosis-prev_est.pdf')
pl.savefig('/homes/peterhm/gbd/book/applications/cirrhosis-prev_est.png')
