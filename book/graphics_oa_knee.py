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
        model = dismod3.data.load('/home/j/Project/dismod/output/dm-38895') 
    except:
        model = dismod3.data.load('/home/j/Project/dismod/dismod_status/prod/dm-38895')
    model.keep(areas=['europe_western'], sexes=['female'])
    model.parameters['m_with'] = {}
    model.input_data = model.input_data.drop(['x_cv_dx_oa_self-report','x_cv_dx_radiographic_only','x_cv_dx_self-reported_pain','x_cv_dx_symptomatic_only'], axis=1)
    
    model.input_data['effective_sample_size'] = 10000
    model.parameters['i']['parameter_age_mesh'] = [0, 30, 31, 35, 40, 65, 100]
    return model

best_model = load_new_model()

# figure oa_knee-knots
output = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_oa_knee.csv')
pl.figure(**book_graphics.full_page_params)

param_list = [dict(type='p', title='(a)', ylabel='Prevalence (%)', yticks=([0, .2, .4, .1, .3], [0, 20, 40, 10, 30]), axis=[25,105,-.045,.45], loc='upper left'),
          dict(type='i', title='(b)', ylabel='Incidence\n(per 1000 PY)', yticks=([0, .003, .006, .009, .012], [0, 3, 6, 9, 12]), axis=[-5,105,-.0015,.015], loc='upper right'),
          dict(type='f', title='(c)', ylabel='Excess mortality\n (per 10,000 PY)', yticks=([0, .0003, .0006, .0009, .0012], [0, 3, 6, 9, 12]), axis=[-5,105,-.00015,.0015], loc='upper left' ),
          #dict(type='r', title='(c)', ylabel='Excess mortality (Per 100 PY)', yticks=([0, .002, .004, .006, .008], [0, 2, 4, 6, 8]), axis=[-5,105,-.0018,.009], loc='upper left' ),
          #dict(type='m_with', title='(d)', ylabel='With-condition mortality \n (per 100 PY)', yticks=([0, .1, .2, .3, .4], [0, 10, 20, 30, 40]), axis=[-5,105,-.045,.45], loc='upper left'),
          ]

for i, params in enumerate(param_list):
    ax = pl.subplot(2,2,i+1)
    dismod3.graphics.plot_data_bars(best_model.get_data(params['type']), color='grey') 

    pl.plot(pl.arange(101), pl.array(output['0k_'+params['type']]), 'k:', linewidth=3, label='{}')
    pl.plot(pl.arange(101), pl.array(output['1k_'+params['type']]), 'k--', linewidth=3, label='{35}')
    pl.plot(pl.arange(101), pl.array(output['2k_'+params['type']]), 'k-', linewidth=3, label='{31, 35}')
    
    pl.xlabel('Age (years)')
    pl.ylabel(params['ylabel']+'\n\n', ha='center')
    pl.yticks(*params.get('yticks', ([0, .025, .05], [0, 2.5, 5])))
    pl.axis(params.get('axis', [-5,105,-.005,.06]))
    subtitle(params['title'])
    pl.grid()
pl.subplots_adjust(hspace=.35)
pl.subplots_adjust(wspace=.35)

pl.legend(bbox_to_anchor=(.42, 0, .27, .47), bbox_transform=pl.gcf().transFigure, fancybox=True, shadow=True, title='Additional knots at:')
pl.savefig('book/graphics/oa_knee-knots.pdf')
pl.savefig('book/graphics/oa_knee-knots.png')

# figure oa_knee-i_prior
pl.figure(**book_graphics.half_page_params)

param_list = [dict(type='p', title='(a)', ylabel='Prevalence (%)', yticks=([0, .12, .24, .36, .48], [0, 12, 24, 36, 48]), axis=[25,105,-.054,.54], loc='upper left'),
          dict(type='i', title='(b)', ylabel='Incidence \n (per 1000 PY)', yticks=([0, .004, .008, .012, .016], [0, 4, 8, 12, 16]), axis=[-5,105,-.0027,.027], loc='upper left'),
          #dict(type='f', title='(c)', ylabel='Excess mortality (Per 100 PY)', yticks=([0, .002, .004, .006, .008], [0, 2, 4, 6, 8]), axis=[-5,105,-.0018,.009], loc='upper left' ),
          #dict(type='r', title='(c)', ylabel='Excess mortality (Per 100 PY)', yticks=([0, .002, .004, .006, .008], [0, 2, 4, 6, 8]), axis=[-5,105,-.0018,.009], loc='upper left' ),
          #dict(type='m_with', title='(d)', ylabel='With-condition mortality (Per 100 PY)', yticks=([0, .1, .2, .3, .4], [0, 10, 20, 30, 40]), axis=[-5,105,-.045,.45], loc='upper left'),
          ]

for i, params in enumerate(param_list):
    ax = pl.subplot(1,2,i+1)
    dismod3.graphics.plot_data_bars(best_model.get_data(params['type']), color='grey') 

    pl.plot(pl.arange(101), pl.array(output['2k_b30_'+params['type']]), 'k--', linewidth=3, label='$i(a)=$ $0$, for $a <$ $30$')
    pl.plot(pl.arange(101), pl.array(output['2k_'+params['type']]), 'k-', linewidth=3, label='$i(a)=$ $0$, for $a <$ $30$\nand $a >$ $99$')
    
    pl.xlabel('Age (years)')
    pl.ylabel(params['ylabel']+'\n\n', ha='center')
    pl.yticks(*params.get('yticks', ([0, .005, .01, .015, .02], [0, 5, 10, 15, 20])))
    pl.axis(params.get('axis', [-5,105,-.0022,.022]))
    subtitle(params['title'])
    pl.grid()
pl.subplots_adjust(wspace=.35, bottom=.15)
pl.legend(bbox_to_anchor=(.42, 0, .5, .94), bbox_transform=pl.gcf().transFigure, fancybox=True, shadow=True)

pl.savefig('book/graphics/oa_knee-i_prior.pdf')
pl.savefig('book/graphics/oa_knee-i_prior.png')

pl.show()