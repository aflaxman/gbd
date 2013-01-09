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
    try:
        model = dismod3.data.load('/home/j/Project/dismod/output/dm-40552')
    except:
        model = dismod3.data.load('/home/j/Project/dismod/dismod_status/prod/dm-40552') # data - 32281, newest model 39661
    model.keep(areas=['europe_western'])
    # delete regional prevalence points
    drop_pts = model.input_data[(model.input_data['area']=='europe_western') & (model.input_data['data_type']=='p')]
    model.input_data = model.input_data.drop(drop_pts.index)
    # fill in csmr covariate nan as zero
    model.input_data['x_ihme_fao_stimulants_kcal_26oct11'] = model.input_data['x_ihme_fao_stimulants_kcal_26oct11'].fillna([0])
    model.input_data['x_smoking_prev'] = model.input_data['x_smoking_prev'].fillna([0])
    return model

best_model = load_new_model()
output = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-parkinsons.csv')

# figure parkinsons-data
pl.figure(**book_graphics.full_page_params)

pl.subplot(2,2,1)
dismod3.graphics.plot_data_bars(best_model.get_data('p'))
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .01, .02], [0, 1, 2])
pl.axis([60,101,-0.001,.025])
subtitle('(a)')
pl.grid()

pl.subplot(2,2,2)
dismod3.graphics.plot_data_bars(best_model.get_data('i'))
pl.xlabel('Age (years)')
pl.ylabel('Incidence \n(per 10,000 PY)\n\n', ha='center')
pl.yticks([0, .001,.002, .003, .004], [0, 1, 2, 3, 4]) 
pl.axis([60,104,-.0003,.0055])
subtitle('(b)')
pl.grid()

pl.subplot(2,2,3)
dismod3.graphics.plot_data_bars(best_model.get_data('csmr'))
pl.xlabel('Age (years)')
pl.ylabel('Cause-specific mortality \n(per 1000 PY)\n\n', ha='center')
pl.yticks([0, .001,.002, .003, .004], [0, 1, 2, 3, 4])
pl.axis([60,104,-.0002,.005])
subtitle('(c)')
pl.grid()

pl.subplot(2,2,4)
dismod3.graphics.plot_data_bars(best_model.get_data('smr'))
pl.xlabel('Age (years)')
pl.ylabel('Standardized \nmortality ratio\n\n', ha='center')
pl.yticks([1, 2, 3,4, ], [1, 2,3, 4])
pl.axis([60,104,.3,4.5])
subtitle('(d)')
pl.subplots_adjust(hspace=.35,wspace=.35)
pl.grid()

pl.savefig('parkinsons-data.pdf')

# parkinsons-best
pl.figure(**book_graphics.full_page_params)

param_list = [dict(type='p', title='(a)', ylabel='Prevalence (%)', yticks=([0, .01, .02], [0, 1, 2]), axis=[60,101,-0.001,.025]),
          dict(type='i', title='(b)', ylabel='Incidence \n(per 1000 PY)', yticks=([0, .001,.002, .003, .004], [0, 1, 2, 3, 4]), axis=[60,104,-.0003,.0055]),
          dict(type='pf', title='(c)', ylabel='Cause-specific mortality \n(per 1000 PY)', yticks=([0, .001,.002, .003, .004], [0, 1, 2, 3, 4]), axis=[60,104,-.0002,.005]),
          dict(type='smr', title='(d)', ylabel='Standardized \nmortality ratio', yticks=([1, 2, 3,4, ], [1, 2,3, 4]), axis=[60,104,.3,4.5]),
          ]

for i, params in enumerate(param_list):
    ax = pl.subplot(2,2,i+1)
    if params['type'] == 'pf': dismod3.graphics.plot_data_bars(best_model.get_data('csmr'), color='grey')
    else: dismod3.graphics.plot_data_bars(best_model.get_data(params['type']), color='grey')
    
    pl.plot(pl.arange(101), pl.array(output['model_'+params['type']]), 'k-', linewidth=2, label='Posterior Mean')
    pl.plot(pl.arange(101), pl.array(output['model_'+params['type']+'l']), 'k-', linewidth=1, label='95% HPD interval')
    pl.plot(pl.arange(101), pl.array(output['model_'+params['type']+'u']), 'k-', linewidth=1)
    
    pl.xlabel('Age (years)')
    pl.ylabel(params['ylabel']+'\n\n', ha='center')
    pl.axis(params.get('axis', [-5,105,-.005,.06]))
    pl.yticks(*params.get('yticks', ([0, .025, .05], [0, 2.5, 5])))
    subtitle(params['title'])
    pl.grid()
    
pl.subplots_adjust(hspace=.35, wspace=.35)

pl.savefig('parkinsons-best.pdf')

pl.show()
