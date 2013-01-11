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
    try:
        model = dismod3.data.load('/home/j/Project/dismod/output/dm-32458')
    except:
        model = dismod3.data.load('/home/j/Project/dismod/dismod_status/prod/dm-32458')
    model.keep(areas=['europe_western'], sexes=['male', 'total'])
    model.input_data = model.input_data.drop(['x_lnASDR_B03.4.3'],1)
    model.input_data.ix[model.get_data('i').index, 'area'] = 'all'
    model.parameters['i']['fixed_effects']['x_sex'] = dict(dist='Normal', mu=0, sigma=.0001)
    return model

best_model = load_new_model()

output = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_af.csv')

# figure af-data
pl.figure(**book_graphics.half_page_params)

pl.subplot(1,2,1)
dismod3.graphics.plot_data_bars(best_model.get_data('p'))
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .05, .1, .15,  .2], [0, 5, 10, 15, 20])
my_axis(.22)
book_graphics.subtitle('(a)')


pl.subplot(1,2,2)
dismod3.graphics.plot_data_bars(best_model.get_data('i'))
pl.xlabel('Age (years)')
pl.ylabel('Incidence (per 1000 PY)')
pl.yticks([0, .001, .002, .003, .004], [0, 1, 2, 3, 4,])
my_axis(.0045)
book_graphics.subtitle('(b)')


pl.subplots_adjust(wspace=.35, bottom=.14)

pl.savefig('book/graphics/af-data.pdf')
pl.savefig('book/graphics/af-data.png')

# figure af-mp_v_hetero_srt_p
pl.figure(**book_graphics.three_quarter_page_params)
x = best_model.parameters['p']['parameter_age_mesh']

pl.subplot(1,2,1)    
dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')
pl.plot(pl.arange(101), pl.array(output['as_sp']), 'k-', linewidth=2, label='Age-standardizing')
pl.plot(x, pl.array(output['as_sp_l'])[x], 'k-', linewidth=1, label='95% HPD interval')
pl.plot(x, pl.array(output['as_sp_u'])[x], 'k-', linewidth=1)

pl.xlabel('Age (Years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .05, .1, .15,  .2], [0, 5, 10, 15, 20])
my_axis(.22)
book_graphics.subtitle('(a)')

pl.legend(loc='upper center', fancybox=True, shadow=True, bbox_to_anchor=(.5,-.13))

pl.subplot(1,2,2)    
dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')
pl.plot(pl.arange(101), pl.array(output['m_sp']), 'k-', linewidth=2, label='Midpoint')
pl.plot(x, pl.array(output['m_sp_l'])[x], 'k-', linewidth=1, label='95% HPD interval')
pl.plot(x, pl.array(output['m_sp_u'])[x], 'k-', linewidth=1)
    
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .05, .1, .15,  .2], [0, 5, 10, 15, 20])
my_axis(.22)
book_graphics.subtitle('(b)')

pl.legend(loc='upper center', fancybox=True, shadow=True, bbox_to_anchor=(.5,-.13))

pl.subplots_adjust(top=.99, bottom=.27, wspace=.35)

pl.savefig('book/graphics/af-mp_v_hetero_srt_p.pdf')
pl.savefig('book/graphics/af-mp_v_hetero_srt_p.png')

# figure af-mp_v_hetero_srt_i
pl.figure(**book_graphics.three_quarter_page_params)

x = best_model.parameters['i']['parameter_age_mesh']
    
pl.subplot(1,2,1)    
dismod3.graphics.plot_data_bars(best_model.get_data('i'), color='grey')
pl.plot(pl.arange(101), pl.array(output['as_si']), 'k-', linewidth=2, label='Age-standardizing')
pl.plot(x, pl.array(output['as_si_l'])[x], 'k-', linewidth=1, label='95% HPD interval')
pl.plot(x, pl.array(output['as_si_u'])[x], 'k-', linewidth=1)
    
pl.xlabel('Age (years)')
pl.ylabel('Incidence (per 1000 PY)')
pl.yticks([0, .001, .002, .003,  .004], [0, 1, 2, 3, 4])  
my_axis(.0045)
book_graphics.subtitle('(a)')

pl.legend(loc='upper center', fancybox=True, shadow=True, bbox_to_anchor=(.5,-.13))

pl.subplot(1,2,2)    
dismod3.graphics.plot_data_bars(best_model.get_data('i'), color='grey')
pl.plot(pl.arange(101), pl.array(output['m_si']), 'k-', linewidth=2, label='Midpoint')
pl.plot(x, pl.array(output['m_si_l'])[x], 'k-', linewidth=1, label='95% HPD interval')
pl.plot(x, pl.array(output['m_si_u'])[x], 'k-', linewidth=1)
    
pl.xlabel('Age (years)')
pl.ylabel('Incidence (per 1000 PY)')
pl.yticks([0, .001, .002, .003,  .004], [0, 1, 2, 3, 4])  
my_axis(.0045)
book_graphics.subtitle('(b)')

pl.legend(loc='upper center', fancybox=True, shadow=True, bbox_to_anchor=(.5,-.13))

pl.subplots_adjust(top=.99, bottom=.27, wspace=.35)

pl.savefig('book/graphics/af-mp_v_hetero_srt_i.pdf')
pl.savefig('book/graphics/af-mp_v_hetero_srt_i.png')

# figure af-best_model
pl.figure(**book_graphics.three_quarter_page_params)

param_list = [dict(type='p', title='(a)', ylabel='Prevalence (%)', yticks=([0, .05, .1, .15,  .2], [0, 5, 10, 15, 20]), axis=[-5,105,-0.022,.22]),
          dict(type='i', title='(b)', ylabel='Incidence (per 1000 PY)', yticks=([0, .001, .002, .003, .004], [0, 1, 2, 3, 4,]), axis=[-5,105,-.00045,.0045]),
          #dict(type='r', title='(c)', ylabel='Remission (Per 100 PY)', yticks=([0, .025, .05], [0, 2.5, 5]), axis=[-5,105,-.005,.05]),
          #dict(type='m_with', title='(d)', ylabel='With-condition mortality (per 100 PY)', yticks=([0, .05, .1]), axis=[-5,105,-.01,.1]),
          ]

for i, params in enumerate(param_list):
    ax = pl.subplot(1,2,i+1)
    dismod3.graphics.plot_data_bars(best_model.get_data(params['type']), color='grey')
    x = best_model.parameters[params['type']]['parameter_age_mesh']
    
    pl.plot(pl.arange(101), pl.array(output['as_c'+params['type']]), 'k-', linewidth=2, label='Posterior mean')
    pl.plot(x, pl.array(output['as_c'+params['type']+'_l'])[x], 'k-', linewidth=1, label='95% HPD interval')
    pl.plot(x, pl.array(output['as_c'+params['type']+'_u'])[x], 'k-', linewidth=1)
    
    pl.xlabel('Age (years)')
    pl.ylabel(params['ylabel'])
    pl.yticks(*params.get('yticks', ([0, .025, .05], [0, 2.5, 5])))
    pl.axis(params.get('axis'))
    book_graphics.subtitle(params['title'])
    
pl.legend(loc='upper center', bbox_to_anchor=(-.2,-.13), fancybox=True, shadow=True)
    
pl.subplots_adjust(top=.99, bottom=.27, wspace=.35)

pl.savefig('book/graphics/af-best_model.pdf')
pl.savefig('book/graphics/af-best_model.png')

# figure af-mp_v_hetero
pl.figure(**book_graphics.half_page_params)

param_list = [dict(type='p', title='(a)', ylabel='Prevalence (%)', yticks=([0, .03, .06, .09,  .12], [0, 3, 6, 9, 12]), axis=[-5,105,-0.015,.15]),
          dict(type='i', title='(b)', ylabel='Incidence (per 1000 PY)', yticks=([0, .001, .002, .003, .004], [0, 1, 2, 3, 4,]), axis=[-5,105,-.00045,.0045]),
          ]

for i, params in enumerate(param_list):
    ax = pl.subplot(1,2,i+1)

    pl.plot(pl.arange(101), pl.array(output['as_c'+params['type']]), 'k-', linewidth=2, label='Age-standardizing')
    pl.plot(pl.arange(101), pl.array(output['m_c'+params['type']]), 'k--', linewidth=2, label='Midpoint')
    
    pl.xlabel('Age (years)')
    pl.ylabel(params['ylabel'])
    pl.yticks(*params.get('yticks', ([0, .025, .05], [0, 2.5, 5])))
    pl.axis(params.get('axis', [-5,105,-.005,.06]))
    book_graphics.subtitle(params['title'])

pl.legend(loc='upper center', bbox_to_anchor=(-.2,-.13), fancybox=True, shadow=True)    
pl.subplots_adjust(top=.99, bottom=.27, wspace=.35)

pl.savefig('book/graphics/af-mp_v_hetero.pdf')
pl.savefig('book/graphics/af-mp_v_hetero.png')

pl.show()