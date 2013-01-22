# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Set up the dismod/python environment

# <codecell>

# fast for development, slow for final draft
fast=False
if fast:
    iter=100
    burn=0
    thin=1
else:
    iter=10000
    burn=5000
    thin=5

# <codecell>

import sys
sys.path += ['../gbd', '../gbd/book', '../dm3-computation_only/', '../dm3-computation_only/book']
sys.path += ['..', '../..']
import pylab as pl
import pymc as mc
import pandas

import dismod3
reload(dismod3)

import book_graphics
reload(book_graphics)

# set font
book_graphics.set_font()


# Graphical summary of data

model = dismod3.data.load('/home/j/Project/dismod/notebooks/models/bipolar_orig')

# <codecell>

def my_axis(ymax):
    axis([-5,105,-ymax/10.,ymax])

model = dismod3.data.load('/home/j/Project/dismod/notebooks/models/bipolar')
df = model.get_data('p')

pl.figure(**book_graphics.half_page_params)
for i in range(2):
    pl.subplot(1,2,1+i)
    dismod3.graphics.plot_data_bars(df[df['x_cv_past_year'] == i])
    pl.xlabel('Age (years)')
    pl.ylabel('Prevalence (%)')
    pl.yticks([0, .01, .02, .03, .04], [0, 1, 2, 3, 4])
    pl.axis([-5,105,-.0045, .045])
    if i == 0: book_graphics.subtitle('(a) Past-year prevalence')
    else: book_graphics.subtitle('(b) Past-month prevalence')
    #title('cv_past_year = %d'%i)

pl.subplots_adjust(wspace=.35, bottom=.14)    
pl.savefig('book/graphics/bipolar-data-by-cv.pdf')


# Compare alternative reference values
# -------------------------------------

model = dismod3.data.load('/home/j/Project/dismod/notebooks/models/bipolar')
#model.keep(areas=['super-region_0'], sexes=['male', 'total'], end_year=1997)

# remove expert prior on pyp effect
model.parameters['p']['fixed_effects'].pop('x_cv_past_year')

model.vars += dismod3.ism.consistent(model)

dismod3.fit.fit_consistent(model, iter=iter, thin=thin, burn=burn)

# <codecell>

py_ref = model

# <codecell>

model = dismod3.data.load('/home/j/Project/dismod/notebooks/models/bipolar')
#model.keep(areas=['super-region_0'], sexes=['male', 'total'], end_year=1997)
model.output_template['x_cv_past_year'] = 1.
# remove expert prior on pyp effect
model.parameters['p']['fixed_effects'].pop('x_cv_past_year')

model.vars += dismod3.ism.consistent(model)

dismod3.fit.fit_consistent(model, iter=iter, thin=thin, burn=burn)

# <codecell>

pm_ref = model

# <codecell>

pl.figure(**book_graphics.three_quarter_page_params)
param_list = [dict(model=py_ref, linestyle='-', label='Past year'),
              dict(model=pm_ref, linestyle='--', label='Past month')]             
for params in param_list:
    model = params['model']
    p = dismod3.covariates.predict_for(model, model.parameters['p'], 'all', 'total', 'all', 'europe_western', 'male', 1990, 1., model.vars['p'], 0., 1.).T
    i = dismod3.covariates.predict_for(model, model.parameters['i'], 'all', 'total', 'all', 'europe_western', 'male', 1990, 1., model.vars['i'], 0., 1.).T
    pl.subplot(1,2,1)
    pl.plot(pl.arange(101), p.mean(axis=1), 'k', linestyle=params['linestyle'], linewidth=2)#, label=params['label'])
    pl.subplot(1,2,2)
    pl.plot(pl.arange(101), i.mean(axis=1), 'k', linestyle=params['linestyle'], linewidth=2)#, label=params['label'])
    
pl.subplot(1,2,1)
pl.xlabel('Age (years)', fontsize='xx-large')
pl.ylabel('Prevalence (%)', fontsize='xx-large')
pl.yticks([0, .005, .01, .015, .02], [0, 0.5, 1.0, 1.5, 2.0], fontsize='x-large')
pl.axis([-5,105,-.0022, .022])
book_graphics.subtitle('(a)')

pl.subplot(1,2,2)
pl.xlabel('Age (years)', fontsize='xx-large')
pl.ylabel('Incidence \n (per 10,000 PY)'+'\n\n', ha='center',fontsize='xx-large')
pl.yticks([0, .0005, .001, .0015, .002], [0, 5, 10, 15, 20], fontsize='x-large')
pl.plot([-10],[-10],'k-', label='Past year')
pl.plot([-10],[-10],'k--', label='Past month')
pl.axis([-5,105,-.00023, .0023])
book_graphics.subtitle('(b)')
pl.subplots_adjust(wspace=.35, bottom=.14)

pl.legend(loc='upper center', fancybox=True, shadow=True, bbox_to_anchor=(-.2,-.2), ncol=2)
pl.subplots_adjust(top=.99, bottom=.27, wspace=.35)

pl.savefig('book/graphics/bipolar-ref-alts.pdf')

pl.show()
