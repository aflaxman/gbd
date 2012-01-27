# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import colorbrewer
colors = ['#%x%x%x' % col for col in colorbrewer.Set1[3]]

import dismod3
reload(dismod3)

model = dismod3.data.ModelData.load(dismod3.settings.JOB_WORKING_DIR % 29820)
model.keep(areas=['super-region_0'])
model.vars += dismod3.ism.age_specific_rate(model, data_type='pf',
                                            reference_area='all', reference_sex='total', reference_year='all')
dismod3.fit.fit_asr(model.vars, iter=30000, burn=15000, thin=15, tune_interval=100) # this should fit all of the asrs available
model.vars.plot_acorr()
model.vars.plot_trace()

dismod3.graphics.plot_one_type(model, model.vars['pf'], {}, 'pf')
pl.show()

