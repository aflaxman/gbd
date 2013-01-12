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
    orig_model = dismod3.data.load('/home/j/Project/dismod/notebooks/models/bipolar_orig') 
    new_model = dismod3.data.load('/home/j/Project/dismod/notebooks/models/bipolar')
    
    nonprev_ix = (orig_model.input_data['data_type'] != 'p')
    df = orig_model.input_data[nonprev_ix]
    
    prev_ix = (new_model.input_data['data_type'] == 'p')
    df = df.append(new_model.input_data[prev_ix], ignore_index=True)
    
    # fill in missing values of x_cv_outlier column
    df['x_cv_outlier'] = df['x_cv_outlier'].fillna(0.)
    
    new_model.input_data = df
    new_model.keep(areas=['north_america_high_income']) 
    #orig_model.parameters['p'] = new_model.parameters['p']
    #model.parameters['i']['level_value'] = dict(value=0, age_before=10, age_after=100)
    return new_model

best_model = load_new_model()

output = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_bipolar.csv')

# figure bipolar-data
pl.figure(**book_graphics.full_page_params)

pl.subplot(2,2,1)
dismod3.graphics.plot_data_bars(best_model.get_data('p'))
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .01, .02, .03, .04], [0, 1, 2, 3, 4])
my_axis(.045)
book_graphics.subtitle('(a)')


pl.subplot(2,2,2)
dismod3.graphics.plot_data_bars(best_model.get_data('i'))
pl.xlabel('Age (years)')
pl.ylabel('Incidence \n (per 10,000 PY)'+'\n\n', ha='center')
pl.yticks([0, .0003, .0006, .0009, .0012], [0, 3, 6, 9, 12])
my_axis(.0014)
book_graphics.subtitle('(b)')


pl.subplot(2,2,3)
dismod3.graphics.plot_data_bars(best_model.get_data('r'))
pl.xlabel('Age (years)')
pl.ylabel('Remission \n (per 100 PY)'+'\n\n', ha='center')
pl.yticks([0, .01, .02, .03, .04], [0, 1, 2, 3, 4])
my_axis(.045)
book_graphics.subtitle('(c)')


pl.subplot(2,2,4)
dismod3.graphics.plot_data_bars(best_model.get_data('smr'))
pl.xlabel('Age (years)')
pl.ylabel('Standardized \n mortality ratio'+'\n\n', ha='center')
pl.yticks([0, 2, 4, 6, 8])
my_axis(9)
book_graphics.subtitle('(d)')


pl.subplots_adjust(hspace=.35)
pl.subplots_adjust(wspace=.35)

pl.savefig('book/graphics/bipolar-data.pdf')
pl.savefig('book/graphics/bipolar-data.png')

# figure bipolar-bounds
pl.figure(**book_graphics.full_page_params)

pl.plot(pl.arange(101), pl.array(output['no_bounds_p']), 'k-', linewidth=2, label='$h_i(a)$ unrestricted')
pl.plot(pl.arange(101), pl.array(output['best_p']), 'k--', linewidth=2, label='$h_i(a) =$ $0$ for $a <$ $10$')

#dismod3.graphics.plot_data_bars(best_model.get_data('p'))
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .005, .01, .015, .02], [0, 0.5, 1.0, 1.5, 2.0])
my_axis(.021)
pl.legend(loc='upper right', fancybox=True, shadow=True)


pl.savefig('book/graphics/bipolar-bounds.pdf')
pl.savefig('book/graphics/bipolar-bounds.png')

# figure bipolar-45_65_100
pl.figure(**book_graphics.full_page_params)

param_list = [dict(age=100, linestyle='-', label='$h_i(a)$ unrestricted'),
              dict(age=65, linestyle='--', label='$h_i(a) =$ $0$ for $a >$ $65$'),
              dict(age=45, linestyle=':', label='$h_i(a) =$ $0$ for $a >$ $45$')]
for params in param_list:
    pl.subplot(2,2,1)
    pl.plot(pl.arange(101), pl.array(output['i'+str(params['age'])+'_p']), 'k', linestyle=params['linestyle'], linewidth=2, label=params['label'])
    pl.subplot(2,2,2)
    pl.plot(pl.arange(101), pl.array(output['i'+str(params['age'])+'_i']), 'k', linestyle=params['linestyle'], linewidth=2, label=params['label'])
    pl.subplot(2,2,3)
    pl.plot(pl.arange(101), pl.array(output['i'+str(params['age'])+'_r']), 'k', linestyle=params['linestyle'], linewidth=2, label=params['label'])
    pl.subplot(2,2,4)
    pl.plot(pl.arange(101), pl.array(output['i'+str(params['age'])+'_f']), 'k', linestyle=params['linestyle'], linewidth=2, label=params['label'])
    
pl.subplot(2,2,1)
#dismod3.graphics.plot_data_bars(model.get_data('p'))
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .005, .01, .015, .02], [0, .5, 1, 1.5, 2])
my_axis(.022)
book_graphics.subtitle('(a)')

pl.subplot(2,2,2)
pl.xlabel('Age (years)')
pl.ylabel('Incidence \n (per 10,000 PY)'+'\n\n', ha='center')
pl.yticks([0, .0005, .001, .0015, .0020], [0, 5, 10, 15, 20])
my_axis(.0022)
book_graphics.subtitle('(b)')

pl.subplot(2,2,3)
pl.xlabel('Age (years)')
pl.ylabel('Remission \n (per 1000 PY)'+'\n\n', ha='center')
pl.yticks([0, .007, .014, .021, .028], [0, 7, 14, 21, 28])
my_axis(.032)
book_graphics.subtitle('(c)')

pl.subplot(2,2,4)
pl.xlabel('Age (years)')
pl.ylabel('Excess mortality \n (per 100 PY)'+'\n\n', ha='center')
pl.yticks([0, .05, .1, .15, .2], [0, 5, 10, 15, 20])
my_axis(.22)
book_graphics.subtitle('(d)')

pl.legend(loc='upper center', bbox_to_anchor=(-.2,-.3), fancybox=True, shadow=True, ncol=3)    
pl.subplots_adjust(top=.99, bottom=.18, wspace=.35, hspace=.3)

pl.savefig('book/graphics/bipolar-45_65_100.pdf')
pl.savefig('book/graphics/bipolar-45_65_100.png')

# figure bipolar-0_5_10
pl.figure(**book_graphics.three_quarter_page_params)

param_list = [dict(ub=.1, linestyle='-', label='$h_r(a) <$ $10$'),
              dict(ub=.05, linestyle='--', label='$h_r(a) <$ $5$'),
              dict(ub=0, linestyle=':', label='$h_r(a) =$ $0$')]
for params in param_list:
    pl.subplot(1,2,1)
    pl.plot(pl.arange(101), pl.array(output['r'+str(params['ub'])+'_r']), 'k', linestyle=params['linestyle'], linewidth=2, label=params['label'])
    
    pl.subplot(1,2,2)
    pl.plot(pl.arange(101), pl.array(output['r'+str(params['ub'])+'_f']), 'k', linestyle=params['linestyle'], linewidth=2, label=params['label'])
    
pl.subplot(1,2,1)
pl.xlabel('Age (years)')
pl.ylabel('Remission (per 100 PY)')
pl.yticks([0, .02, .04, .06, .08], [0, 2, 4, 6, 8])
my_axis(.09)
book_graphics.subtitle('(a)')

pl.subplot(1,2,2)
pl.xlabel('Age (years)')
pl.ylabel('Excess mortality \n (per 100 PY)'+'\n\n', ha='center')
pl.yticks([0, .04, .08, .12, .16], [0, 4, 8, 12, 16])
my_axis(.18)
book_graphics.subtitle('(b)')

pl.legend(loc='upper center', bbox_to_anchor=(-.2,-.2), fancybox=True, shadow=True, ncol=3)    
pl.subplots_adjust(top=.99, bottom=.27, wspace=.35)

pl.savefig('book/graphics/bipolar-0_5_10.pdf')
pl.savefig('book/graphics/bipolar-0_5_10.png')

pl.show()