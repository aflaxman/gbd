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
    
def load_new_model():
    id = 39605
    try:
        model = dismod3.data.load('/home/j/Project/dismod/output/dm-39605') 
    except:
        model = dismod3.data.load('/home/j/Project/dismod/dismod_status/prod/dm-39605')
    model.keep(areas=['asia_central'], sexes=['male', 'total'])
    # add sex FE
    model.parameters['f']['fixed_effects']['x_sex'] = dict(dist='Normal', mu=0., sigma=.0001)
    model.parameters['p']['fixed_effects']['x_sex'] = dict(dist='Normal', mu=0., sigma=.0001)
    # drop covartiates
    model.input_data = model.input_data.drop(['x_alcohol_lpc','x_lnASDR_B08.2.2','x_wmhs'], axis=1)
    # make excess mortality 0 from ages [0, 10]
    #model.parameters['f']['level_value'] = {'age_after': 100, 'age_before': 10, 'value': 0.0}
    model.parameters['f']['parameter_age_mesh'] = [0, 100]
    return model

csmr_model = load_new_model()

# figure alcohol-data
pl.figure(**book_graphics.full_page_params)

pl.subplot(2,2,1)
dismod3.graphics.plot_data_bars(csmr_model.get_data('p'))
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .03, .06, .09, .12], [0, 3, 6, 9, 12])
my_axis(.13)
book_graphics.subtitle('(a)')


pl.subplot(2,2,2)
dismod3.graphics.plot_data_bars(csmr_model.get_data('i'))
pl.xlabel('Age (years)')
pl.ylabel('Incidence (per PY)')
pl.yticks([0, 3, 6, 9, 12])
my_axis(13)
book_graphics.subtitle('(b)')


pl.subplot(2,2,3)
dismod3.graphics.plot_data_bars(csmr_model.get_data('csmr'))
pl.xlabel('Age (years)')
pl.ylabel('Cause-specific mortality \n (per 100,000 PY)'+'\n\n', ha='center')
pl.yticks([0, .00005, .0001, .00015, .00020], [0, 5, 10, 15, 20])
my_axis(.00021)
book_graphics.subtitle('(c)')


pl.subplot(2,2,4)
dismod3.graphics.plot_data_bars(csmr_model.get_data('f'))
pl.xlabel('Age (years)')
pl.ylabel('Excess mortality \n (per 1000 PY)'+'\n\n', ha='center')
pl.yticks([0, .013, .026, .039, .052], [0, 13, 26, 39, 52])
my_axis(.055)
book_graphics.subtitle('(d)')


pl.subplots_adjust(hspace=.35)
pl.subplots_adjust(wspace=.35)

pl.savefig('book/graphics/alcohol-data.pdf')
pl.savefig('book/graphics/alcohol-data.png')    

# figure alcohol-overlay
pl.figure(**book_graphics.full_page_params)
pred = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-alcohol.csv')

param_list = [dict(type='p', title='(a)', ylabel='Prevalence (%)', yticks=([0, .01, .02, .03, .04], [0, 1, 2, 3, 4]), axis=[-5,105,-0.0045,.045]),
          #dict(type='i', title='(b)', ylabel='Incidence (Per 100 PY)', yticks=([0, .03, .06, .09, .12], [0, 3, 6, 9, 12]), axis=[-5,105,-.014,.14]),
          dict(type='pf', title='(b)', ylabel='Cause-specific \n mortality \n (per 100,000 PY)', yticks=([0, .0002, .0004, .0006, .0008], [0, 2, 4, 6, 8]), axis=[-5,105,-.00009,.0009]),
          dict(type='f', title='(c)', ylabel='Excess mortality \n (per 100 PY) \n', yticks=([0, .01, .02, .03, .04], [0, 1, 2, 3, 4]), axis=[-5,105,-.0045,.045]),
          #dict(type='rr', title='(e)', ylabel='Relative Risk (Per 1)', yticks=([0, 1, 10], [0, 1, 10]), axis=[-5,105,-.1,10.1]),
        ]

for i, params in enumerate(param_list):
    ax = pl.subplot(2,2,i+1)

    pl.plot(pl.arange(101), pl.array(pred['csmr_'+params['type']]), 'k-', linewidth=3, label='Posterior Mean, $h_{f^{ \prime \prime}}$ unrestricted')
    pl.plot(pl.arange(101), pl.array(pred['pf_'+params['type']]), 'k--', linewidth=3, label='Posterior Mean, $h_{f^{ \prime \prime}} = 0$')
    
    pl.xlabel('Age (years)')
    if params['type']=='pf': pl.ylabel(params['ylabel']+'\n\n\n', ha='center')
    else: pl.ylabel(params['ylabel']+'\n\n', ha='center')
    pl.axis(params.get('axis', [-5,105,-.005,.06]))
    book_graphics.subtitle(params['title'])
    
    pl.yticks(*params.get('yticks', ([0, .025, .05], [0, 2.5, 5])))
    if i ==2: pl.legend(loc='upper right', bbox_to_anchor=(2.34,1.03), fancybox=True, shadow=True) 
    
pl.subplots_adjust(hspace=.35)
pl.subplots_adjust(wspace=.45)

pl.savefig('book/graphics/alcohol-overlay.pdf')
pl.savefig('book/graphics/alcohol-overlay.png')

pl.show()