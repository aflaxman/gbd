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
        model = dismod3.data.load('/home/j/Project/dismod/output/dm-39098')
    except:
        model = dismod3.data.load('/home/j/Project/dismod/dismod_status/prod/dm-39098')
    # remove covariates
    model.input_data = model.input_data.drop(['x_LDI_id_Updated_7July2011', 'x_ihme_health_system_access_19jul2011'], 1)
    
    # remove increasing age pattern, to make data speak for itself
    #model.parameters['i']['increasing'] = {'age_end': 0, 'age_start': 0}
    #model.parameters['i']['smoothness']['amount'] = 'Moderately'
    return model    

output = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_pancreatitis.csv')    
# figure pancreatitis-we_data
we = load_new_model()
we.keep(areas=['europe_western'], sexes=['male'])

pl.figure(**book_graphics.full_page_params)

dismod3.graphics.plot_data_bars(we.get_data('i'), color='grey')

x=we.parameters['i']['parameter_age_mesh']
pl.plot(pl.arange(101), pl.array(output['we']), 'k-', linewidth=3, label='Posterior')
pl.plot(x, pl.array(output['we_l'])[x], 'k-', linewidth=1, label='95% HPD interval')
pl.plot(x, pl.array(output['we_u'])[x], 'k-', linewidth=1)
        
pl.xlabel('Age (years)')
pl.ylabel('Incidence (per 1000 PY)') #+'\n\n', ha='center')
my_axis(.0052)
pl.grid()
pl.yticks([0, .001, .002, .003, .004], [0, 1, 2, 3, 4])
pl.legend(loc='upper right', fancybox=True, shadow=True)

pl.savefig('graphics/pancreatitis-we_data.pdf')
pl.savefig('graphics/pancreatitis-we_data.png')

# figure pancreatitis-we_compare
fin_wp = load_new_model()
fin_wp.keep(areas=['FIN'], sexes=['male'])
gbr_wp = load_new_model()
gbr_wp.keep(areas=['CYP'], sexes=['male'])
ndl_wp = load_new_model()
ndl_wp.keep(areas=['NLD'], sexes=['male'])
deu_wp = load_new_model()
deu_wp.keep(areas=['DEU'], sexes=['male'])

model_list = [dict(model=fin_wp, subtitle='(a)', cty='FIN', prior=output['FIN_pr']),
              dict(model=ndl_wp, subtitle='(b)', cty='NLD', prior=output['NLD_pr']),
              dict(model=gbr_wp, subtitle='(c)', cty='CYP', prior=output['CYP_pr']),
              dict(model=deu_wp, subtitle='(d)', cty='DEU', prior=output['DEU_pr'])
              ]
              
pl.figure(**book_graphics.full_page_params)

for i, params in enumerate(model_list):
    model = params['model']
    pl.subplot(2,2,i+1) 
    dismod3.graphics.plot_data_bars(model.get_data('i'), color='grey')
    
    x = model.parameters['i']['parameter_age_mesh']
    pl.plot(pl.arange(101), pl.array(output[params['cty']+'_pr']), 'k--', linewidth=3, label='Empirical prior')
    pl.plot(pl.arange(101), pl.array(output[params['cty']]), 'k-', linewidth=3, label='Posterior')
    pl.plot(x, pl.array(output[params['cty']+'_l'])[x], 'k-', linewidth=1, label='95% HPD interval')
    pl.plot(x, pl.array(output[params['cty']+'_u'])[x], 'k-', linewidth=1)

    pl.xlabel('Age (years)')
    pl.ylabel('Incidence \n (per 1000 PY)'+'\n\n', ha='center')
    my_axis(.005)
    subtitle(params['subtitle'])
    pl.grid()
    pl.yticks([0, .001, .002, .003, .004], [0, 1, 2, 3, 4])
    if params['subtitle'] == '(d)': pl.legend(bbox_to_anchor=(.42, 0, .5, .5), bbox_transform=pl.gcf().transFigure, fancybox=True, shadow=True)
  
pl.subplots_adjust(hspace=.35)
pl.subplots_adjust(wspace=.35)
pl.savefig('graphics/pancreatitis-we_compare.pdf')
pl.savefig('graphics/pancreatitis-we_compare_.png')    

pl.show()