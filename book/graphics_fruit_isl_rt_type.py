import sys
sys.path += ['../gbd', '../gbd/book', '../dm3-computation_only/', '../dm3-computation_only/book']
import pylab as pl
import pymc as mc
import pandas

import dismod3
reload(dismod3)

import book_graphics
reload(book_graphics)

import data_model
reload(data_model)

# set font
book_graphics.set_font()

# axis def
def my_axis(ymax):
    pl.axis([-5,105,-ymax/10.,ymax])


def load_we_new_model():
    # example of predicting out-of-sample with a ln_ASDR covariate
    try:
        model = dismod3.data.load('/home/j/Project/dismod/output/dm-35020')
    except:
        model = dismod3.data.load('/home/j/Project/dismod/dismod_status/prod/dm-35020')
    model.keep(areas=['europe_western'], sexes=['male'], start_year=1997)
    model.input_data = model.input_data.drop(['z_cv_natl_rep','x_cv_diet_assess_method','x_cv_met_suboptimal','x_cv_natl_rep','x_fao_factor1','x_fao_factor2','x_fao_factor4','x_ln_LDI_pc','x_ln_fruits'], 1)
    return model

we = load_we_new_model()

# load data to plot
age_pred = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_age_pred.csv', index_col=0)
ui_pred = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_ui_pred.csv', index_col=0)
isl = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_isl.csv', index_col=0)
grc = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_grc.csv', index_col=0)
we = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_we.csv', index_col=0)

GRC_model = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-fruit_GRC_model.csv')
GRC_log_model = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-fruit_GRC_log_model.csv')
GRC_norm_model = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-fruit_GRC_norm_model.csv')

ISL_model = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-fruit_ISL_model.csv')
ISL_log_model = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-fruit_ISL_log_model.csv')
ISL_norm_model = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-fruit_ISL_norm_model.csv')

GRC_data = pandas.DataFrame(pl.zeros((len(GRC_model.columns),3)), columns=['0','1','2'])
ISL_data = pandas.DataFrame(pl.zeros((len(ISL_model.columns),3)), columns=['0','1','2'])

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

for m,model in enumerate([GRC_model,GRC_log_model,GRC_norm_model]):
    for i in list(model.columns):
        GRC_data.ix[int(i),str(m)] = pl.dot(model[i], age_weights[:,1])
for m,model in enumerate([ISL_model,ISL_log_model,ISL_norm_model]):
    for i in list(model.columns):
        ISL_data.ix[int(i),str(m)] = pl.dot(model[i], age_weights[:,1])

# plot 1
labeling = {'GRC':{'x':10,'y':.05,'text':'Greece'}, 
            'ISL':{'x':80,'y':.005,'text':'Iceland'}}
            
pl.figure(**book_graphics.full_page_params)
knots = list(ui_pred.index)

pl.plot(dismod3.graphics.plot_data_bars(we, color='grey', label='Western Europe'))
pl.plot(dismod3.graphics.plot_data_bars(isl, color='black', label='Iceland'))
pl.plot(dismod3.graphics.plot_data_bars(grc, color='black', label='Greece'))
for i in ['GRC', 'ISL']: #in ['we_model', 'we_log_model', 'we_norm_model']:
    pred = pl.array(age_pred['we_model_' + i])
    we_hpd_l = pl.array(ui_pred['we_model_' + i + '_l'])
    we_hpd_u = pl.array(ui_pred['we_model_' + i + '_u'])
    # plot 1
    pl.plot(pred, 'k-', linewidth=3)#, label=labeling[i]['text'])
    pl.plot(knots, we_hpd_l, 'k-', linewidth=1)#, label=labeling['we_model'][i]['text'])
    pl.plot(knots, we_hpd_u, 'k-', linewidth=1)  
    pl.text(labeling[i]['x'], labeling[i]['y'], labeling[i]['text'], size=16, ha='left', va='top')
pl.xlabel('Age (years)')
pl.ylabel('Consumption (kg/d)')
pl.yticks([0, .015, .03, .045, .06], [0, 0.15, 0.30, 0.45, 0.6])
my_axis(.075)

pl.savefig('book/graphics/fruit-we_rate_type.pdf')
pl.savefig('book/graphics/fruit-we_rate_type.png')

# qq plot distribution comparison
pl.figure(**book_graphics.full_plus_page_params)
ix = pl.arange(GRC_data.index[0], GRC_data.index[-1], int(GRC_data.index[-1]*.025))

for i,c in enumerate(['GRC', 'ISL']): #in ['we_model', 'we_log_model', 'we_norm_model']:
    pl.subplot(2,2,i+1)
    pl.plot(pl.sort(pl.array(GRC_data['0'])[ix]), pl.sort(pl.array(GRC_data[str(i+1)])[ix]), 'ko')
    pl.plot([-1,1],[-1,1],'k-')
    pl.yticks([.025, .03, .035, .035], [.25, .3, .35, .35])
    pl.xticks([.025, .03, .035, .035], [.25, .3, .35, .35])
    pl.axis([.023, .042, .023, .041])
    pl.xlabel('Negative-binomial')
    if i + 1 == 1: pl.ylabel('Greece logormal')
    elif i + 1 == 2: pl.ylabel('Greece normal')
    
    pl.subplot(2,2,i+3)
    pl.plot(pl.sort(pl.array(ISL_data['0'])[ix]), pl.sort(pl.array(ISL_data[str(i+1)])[ix]), 'ko')
    pl.plot([-1,1],[-1,1],'k-')
    pl.yticks([.006, .008, .010, .012], [.06, .08, .10, .12])
    pl.xticks([.006, .008, .010, .012], [.06, .08, .10, .12])
    pl.axis([.006, .011, .005, .0135])
    pl.xlabel('Negative-binomial')
    if i + 1 == 1: pl.ylabel('Iceland logormal')
    elif i + 1 == 2: pl.ylabel('Iceland normal')

pl.subplots_adjust(top=.99, bottom=.14, wspace=.35, hspace=.25)
    
pl.savefig('book/graphics/fruit-isl_rate_type.pdf')
pl.savefig('book/graphics/fruit-isl_rate_type.png')

pl.show()