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
    
def load_new_model():
    try:
        model = dismod3.data.load('/home/j/Project/dismod/output/dm-34944')
    except:
        model = dismod3.data.load('/home/j/Project/dismod/dismod_status/prod/dm-34944')
    model.keep(areas=['australasia'], sexes=['female'], start_year=2000)
    # seems to be some trouble with missing values in the mx covariates
    model.input_data = model.input_data.drop(['x_mx_conflict_bin', 'x_mx_shock_10_years_bin'], axis=1)
    model.parameters['p']['fixed_effects']['x_sex'] = dict(dist='Constant', mu=0)
    model.parameters['p']['fixed_effects']['x_cv_below7disorders'] = dict(dist='Constant', mu=0)
    return model

best_model = load_new_model()
pm_model = load_new_model()
pm_model.input_data = pm_model.input_data[pm_model.input_data['x_cv_past_year'] == 0]

# figure anxiety-data_by_cv
pl.figure(**book_graphics.half_page_params)

df = best_model.get_data('p')

pl.figure(**book_graphics.half_page_params)

for i in range(2):
    pl.subplot(1,2,1+i)
    dismod3.graphics.plot_data_bars(df[df['x_cv_past_year'] == i])
    pl.xlabel('Age (years)')
    pl.ylabel('Prevalence (%)')
    pl.yticks([0, .07, .14, .21, .28], [0, 7, 14, 21, 28])
    my_axis(.30)
    pl.grid()
    
    if i == 0: subtitle('(a)')
    if i == 1: subtitle('(b)')
    
pl.subplots_adjust(wspace=.35, bottom=.14)
pl.savefig('book/graphics/anxiety-data_by_cv.pdf')
pl.savefig('book/graphics/anxiety-data_by_cv.png')

# figure anxiety-FE
data = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_anxiety.csv')
pl.figure(**book_graphics.full_page_params)

my_plot_data_bars(best_model.get_data('p'), color='grey', label='Period prevalence')
my_plot_data_bars(pm_model.get_data('p'), color='black',label='Point prevalence')

pl.plot(pl.arange(101), pl.array(data['best']),  'k-', linewidth=2, label='All data, with fixed effects')
pl.plot(pl.arange(101), pl.array(data['no_FE']), 'k--', linewidth=2, label='All data, without fixed effects')
pl.plot(pl.arange(101), pl.array(data['pt_p']),  'k:', linewidth=2, label='Point prevalence')

pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .07, .14, .21, .28], [0, 7, 14, 21, 28])
my_axis(.30)
pl.grid()
pl.legend(loc='upper right', fancybox=True, shadow=True)

pl.savefig('book/graphics/anxiety-FE.pdf')
pl.savefig('book/graphics/anxiety-FE.png')

pl.show()