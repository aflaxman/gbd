import sys
sys.path += ['../gbd', '../gbd/book', '../dm3-computation_only/', '../dm3-computation_only/book']
import pylab as pl
import pymc as mc
import pandas
import random

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


best_model = load_USA_new_model()

output = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_USA_we.csv')

USA_model = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-fruit_USA_model.csv')
USA_log_model = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-fruit_USA_log_model.csv')
USA_norm_model = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-fruit_USA_norm_model.csv')

USA_data = pandas.DataFrame(pl.zeros((len(USA_model.columns),3)), columns=['0','1','2'])

# age standardizing
import scikits.statsmodels.lib.io as pd
aw_file = pandas.DataFrame(pd.genfromdta('/home/j/Project/COD/envelope/data/age_weights_final.dta'))
for a in range(len(aw_file['age'])): aw_file['age'][a] = round(aw_file['age'][a],3)
aw_file = pandas.DataFrame(aw_file['weight'], index=list(aw_file['age']))

age_weights = pl.vstack((pl.arange(101), pl.zeros(101))).transpose()
for a in range(101):
    if a == 0: age_weights[a,1] = aw_file.ix[0.0] + aw_file.ix[0.01] + aw_file.ix[0.1]
    elif (a>=1) & (a<5): age_weights[a,1] = aw_file.ix[1.0]/4
    elif (a>=5) & (a<10): age_weights[a,1] = aw_file.ix[5.0]/5
    elif (a>=10) & (a<15): age_weights[a,1] = aw_file.ix[10.0]/5
    elif (a>=15) & (a<20): age_weights[a,1] = aw_file.ix[15.0]/5
    elif (a>=20) & (a<25): age_weights[a,1] = aw_file.ix[20.0]/5
    elif (a>=25) & (a<30): age_weights[a,1] = aw_file.ix[25.0]/5
    elif (a>=30) & (a<35): age_weights[a,1] = aw_file.ix[30.0]/5
    elif (a>=35) & (a<40): age_weights[a,1] = aw_file.ix[35.0]/5
    elif (a>=40) & (a<45): age_weights[a,1] = aw_file.ix[40.0]/5
    elif (a>=45) & (a<50): age_weights[a,1] = aw_file.ix[45.0]/5
    elif (a>=50) & (a<55): age_weights[a,1] = aw_file.ix[50.0]/5
    elif (a>=55) & (a<60): age_weights[a,1] = aw_file.ix[55.0]/5
    elif (a>=60) & (a<65): age_weights[a,1] = aw_file.ix[60.0]/5
    elif (a>=65) & (a<70): age_weights[a,1] = aw_file.ix[65.0]/5
    elif (a>=70) & (a<75): age_weights[a,1] = aw_file.ix[70.0]/5
    elif (a>=75) & (a<80): age_weights[a,1] = aw_file.ix[75.0]/5
    elif (a>=80) & (a<101): age_weights[a,1] = aw_file.ix[80.0]/5

for m,model in enumerate([USA_model,USA_log_model,USA_norm_model]):
    for i in list(model.columns):
        USA_data.ix[int(i),str(m)] = pl.dot(model[i], age_weights[:,1])

# figure 
fig = pl.figure(**book_graphics.three_quarter_page_params)
ax1 = fig.add_subplot(1,2,1)
ax1.plot(dismod3.graphics.plot_data_bars(best_model.get_data('r'), color='grey'))
   
knots = best_model.parameters['r']['parameter_age_mesh']

ax1.plot(pl.array(output['USA_model']), 'k-', linewidth=3, label='Negative-binomial posterior mean')
ax1.plot(knots, pl.array(output['USA_model_l'])[knots], 'k-', linewidth=1, label='Negative-binomial 95% HPD interval')
ax1.plot(knots, pl.array(output['USA_model_u'])[knots], 'k-', linewidth=1)

pl.xlabel('Age (years)')
pl.ylabel('Consumption (kg/d)')
pl.yticks([0, .006, .012, .018, .024], [0, 0.06, 0.12, 0.18, 0.24])
my_axis(.026)
book_graphics.subtitle('(a)')
#pl.legend(loc='upper center', bbox_to_anchor=(.2,-.23), fancybox=True, shadow=True) 

# subset to plot
ix = pl.arange(USA_data.index[0], USA_data.index[-1], int(USA_data.index[-1]*.025))

ax2 = fig.add_subplot(2,4,3)
ax2.plot(pl.sort(pl.array(USA_data['0'])[ix]), pl.sort(pl.array(USA_data['1'])[ix]), 'ko', label='Lognormal')
ax2.plot([-1,1],[-1,1],'k-')
pl.yticks([.006, .008, .010], [.06, .08, .10])
pl.xticks([.006, .008, .010], [.06, .08, .10])
pl.axis([.0055, .0105, .0055, .0105])
pl.ylabel('Lognormal')
book_graphics.subtitle('(b)')

ax3 = fig.add_subplot(2,4,7, sharex=ax2)
ax3.plot(pl.sort(pl.array(USA_data['0'])[ix]), pl.sort(pl.array(USA_data['2'])[ix]), 'ko', label='Normal')
ax3.plot([-1,1],[-1,1],'k-')
pl.yticks([.006, .008, .010], [.06, .08, .10])
pl.xticks([.006, .008, .010], [.06, .08, .10])
pl.axis([.0055, .0105, .0055, .0105])
pl.ylabel('Normal')
book_graphics.subtitle('(c)')

pl.xlabel('Negative-binomial')

pl.subplots_adjust(top=.99, bottom=.27, wspace=.7, hspace=.1, left=.2, right=1.1) 
pl.setp(ax2.get_xticklabels(), visible=False)
#pl.legend(loc='upper center', bbox_to_anchor=(-.2,-.33), fancybox=True, shadow=True, ncol=3)
 
pl.savefig('book/graphics/fruit-rate_type.pdf')
pl.savefig('book/graphics/fruit-rate_type.png')

