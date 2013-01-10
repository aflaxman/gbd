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

# <markdowncell>

# Bipolar disorder
# ================
# 
# The disease modelling of bipolar disorder was based on literature describing it as a chronic illness
# with little or no complete remission {American Psychiatric Association, 2000 #31}. This approach differed
# to the modelling of chronic episodic mood disorders like major depressive disorder where an estimate of
# disease duration rather than remission is required.
# 
# Thirty-two studies were identified for prevalence, covering 22 countries in 11 GBD world regions (Table 1).
# Two studies were identified for incidence, both from the United States of America. Seven studies were
# identified for excess mortality, from 5 high income countries. We found no studies reporting
# on complete remission (equivalent to cure rather than a temporary reduction in symptom levels 
# as clinicians tend to define ‘remission’) from bipolar disorder as defined by GBD. 
# 
# There was considerable variability in the data. For instance, estimates were reported across multiple age
# groups (specific or broad), they were based on different coverage areas (community or national) and they
# were either sex specific or for both sex combined. This data has been summarised in greater detail
# elsewhere {Ferrari, 2011 #1574; Ferrari, 2011 #1581}. Since there was an overall lack of estimates
# across all parameters and they were characterised by considerable variability, certain assumptions were 
# made to guide the modelling process. These were made in consultation with disease experts.  

# <markdowncell>

# Graphical summary of data
# -------------------------

# <codecell>

model = dismod3.data.load('/home/j/Project/dismod/notebooks/models/bipolar_orig')

# <codecell>

def my_axis(ymax):
    axis([-5,105,-ymax/10.,ymax])

def subtitle(s):
    """ title where the panel names appear within each panel"""
    l,r,b,t=pl.axis()
    x = l + (r-l)*.05
    y = t - (t-b)*.05
    pl.text(x, y, s, ha='left', va='top')    
    
# <codecell>

# # # figure(**book_graphics.full_page_params)

# # # subplot(2,2,1)
# # # dismod3.graphics.plot_data_bars(model.get_data('p'))
# # # title('(a)')
# # # xlabel('Age (years)')
# # # ylabel('Prevalence (%)')
# # # yticks([0, .025, .05], [0, 2.5, 5])
# # # my_axis(.06)

# # # subplot(2,2,2)
# # # dismod3.graphics.plot_data_bars(model.get_data('i'))
# # # title('(b)')
# # # xlabel('Age (years)')
# # # ylabel('Incidence (per 1000 PY)')
# # # yticks([0, .0005, .001], [0, .5, 1])
# # # my_axis(.0014)

# # # subplot(2,2,3)
# # # dismod3.graphics.plot_data_bars(model.get_data('r'))
# # # title('(c)')
# # # xlabel('Age (years)')
# # # ylabel('Remission (per 100 PY)')
# # # yticks([0, .025, .05], [0, 2.5, 5])
# # # my_axis(.06)

# # # subplot(2,2,4)
# # # dismod3.graphics.plot_data_bars(model.get_data('smr'))
# # # title('(d)')
# # # xlabel('Age (years)')
# # # ylabel('Standardized mortality ratio')
# # # yticks([0, 5, 10])
# # # my_axis(15)
# # # subplots_adjust(hspace=.3)

# # # savefig('bipolar-data.pdf')

# # # # <codecell>

# # # model = dismod3.data.load('/snfs1/Project/dismod/notebooks/models/bipolar')

# # # # <codecell>

# # # figure(**book_graphics.full_page_params)

# # # ax = subplot(2,2,1)
# # # dismod3.graphics.plot_data_bars(model.get_data('p'))
# # # title('(a)')
# # # xlabel('Age (years)')
# # # ylabel('Prevalence (%)')
# # # yticks([0, .025, .05], [0, 2.5, 5])

# # # subplot(2,2,2,sharex=ax,sharey=ax)
# # # dismod3.graphics.plot_data_bars(model.get_data('i'))
# # # title('(b)')
# # # xlabel('Age (years)')
# # # ylabel('Incidence (per 100 PY)')
# # # yticks([0, .025, .05], [0, 2.5, 5])

# # # subplot(2,2,3,sharex=ax,sharey=ax)
# # # dismod3.graphics.plot_data_bars(model.get_data('r'))
# # # title('(c)')
# # # xlabel('Age (years)')
# # # ylabel('Remission (per 100 PY)')
# # # yticks([0, .025, .05], [0, 2.5, 5])

# # # subplot(2,2,4,sharex=ax,sharey=ax)
# # # dismod3.graphics.plot_data_bars(model.get_data('f'))
# # # title('(d)')
# # # xlabel('Age (years)')
# # # ylabel('Excess mortality (per 100 PY)')
# # # yticks([0, .025, .05], [0, 2.5, 5])
# # # axis([-5,105,-.005,.06])
# # # subplots_adjust(hspace=.3)

# # # #savefig('bipolar-data.pdf')

# # # # <markdowncell>

# # # # Age of onset 
# # # # ------------
# # # # 
# # # # While there is evidence to suggest that bipolar disorder commonly starts in the mid-teens or early twenties, 
# # # # there is still disagreement over a minimum age of onset {Goodwin, 2008 #208}. Even though symptoms can be 
# # # # tracked back to childhood, setting a threshold for diagnosis is difficult given that current diagnostic 
# # # # criteria are based on adult presentation of the disorder.  The literature and expert advice on this issue 
# # # # suggests that although prepubertal bipolar disorder is rare, a distinct phenotype can exist. In order to
# # # # capture as much of the burden attributable to this disorder, a minimum age of onset of 10 years was set for
# # # # prevalence and incidence (figures 1 and 2).

# # # # <codecell>

# # # model = dismod3.data.load('/snfs1/Project/dismod/notebooks/models/bipolar_orig')
# # # #model.keep(areas=['europe_western'], sexes=['male', 'total'], end_year=1997)

# # # # <codecell>

# # # model.vars = dismod3.ism.consistent(model)

# # # # <codecell>

# # # dismod3.fit.fit_consistent(model, iter=iter, thin=thin, burn=burn)

# # # # <codecell>

# # # best = model

# # # # <codecell>

# # # figure(**book_graphics.full_page_params)

# # # param_list = [dict(type='p', title='(a)', ylabel='Prevalence (%)', yticks=([0, .005, .01], [0, .5, 1]), axis=[-5,105,0,.015]),
          # # # dict(type='i', title='(b)', ylabel='Incidence (per 10,000 PY)', yticks=([0, .0004, .0008], [0, 4, 8]), axis=[-5,105,0,.001]),
          # # # dict(type='r', title='(c)', ylabel='Remission (per 100 PY)'),
          # # # dict(type='f', title='(d)', ylabel='Excess mortality (per 100 PY)', axis=[-5,105,0,.1]),
          # # # ]


# # # for i, params in enumerate(param_list):
    # # # ax = subplot(2,2,i+1)
    # # # dismod3.graphics.plot_data_bars(model.get_data(params['type']))
    
    # # # model.pred = dismod3.covariates.predict_for(model, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                                # # # 1., model.vars[params['type']], 0., 1.).T
    
    # # # plot(arange(101), model.pred.mean(axis=1), 'k-', linewidth=2, label='Posterior Mean')
    # # # hpd = mc.utils.hpd(model.pred.T, .05)
    # # # plot(arange(101), hpd[:,0], 'k-', linewidth=1, label='95% HPD interval')
    # # # plot(arange(101), hpd[:,1], 'k-', linewidth=1)

    # # # title(params['title'])
    # # # xlabel('Age (years)')
    # # # ylabel(params['ylabel'])
    # # # yticks(*params.get('yticks', ([0, .025, .05], [0, 2.5, 5])))
    # # # axis(params.get('axis', [-5,105,-.005,.06]))
# # # subplots_adjust(hspace=.35)

# # # savefig('bipolar-zero_before_ten.pdf')

# # # # <markdowncell>

# # # # What happens without the age restrictions on incidence and prevalence?
# # # # ---------------------------------------------------

# # # # <codecell>

# # # model = dismod3.data.load('/snfs1/Project/dismod/notebooks/models/bipolar')
# # # #model.keep(areas=['europe_western'], sexes=['male', 'total'], end_year=1997)

# # # model.parameters['i']['level_value'] = dict(value=0, age_before=0, age_after=100)
# # # model.parameters['p']['level_value'] = dict(value=0, age_before=0, age_after=100)

# # # model.vars += dismod3.ism.consistent(model)

# # # dismod3.fit.fit_consistent(model, iter=iter, thin=thin, burn=burn)

# # # wo_age = model

# # # # <codecell>

# # # no_bounds = dismod3.covariates.predict_for(wo_age, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                          # # # 1., model.vars['p'], 0., 1.).T
# # # bounds = dismod3.covariates.predict_for(best, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                       # # # 1., best.vars['p'], 0., 1.).T

# # # plot(arange(101), no_bounds.mean(axis=1), 'k-', linewidth=2, label='No prior on age of onset p, i')
# # # plot(arange(101), bounds.mean(axis=1), 'k--', linewidth=2, label='Age of onset > 10')
# # # dismod3.graphics.plot_data_bars(model.get_data('p'))
# # # xlabel('Age (years)')
# # # ylabel('Prevalence (%)')
# # # yticks([0, .005, .01], [0, .5, 1])
# # # legend(loc='upper right', fancybox=True, shadow=True)
# # # #axis(params.get('axis', [-5,105,-.005,.06]))

# # # # <markdowncell>

# # # # What about the restriction of age of onset before a certain age?
# # # # ----------------------------------------------------------------
# # # # 
# # # # An upper age restriction of 65 years was also placed on incidence 
# # # # as it led to the most plausible fit to the data after trialling limits of 45, 65, and 100 years. 

# # # # <codecell>

# # # model.parameters['i']

# # # # <codecell>

# # # # TODO: decide if using global data is more/less instructive than europe west/male/1990

# # # bounds = {}

# # # for upper_bound in [45, 65, 100]:
    # # # model = dismod3.data.load('/snfs1/Project/dismod/notebooks/models/bipolar')
    # # # #model.keep(areas=['europe_western'], sexes=['male', 'total'], end_year=1997)

    # # # model.parameters['i']['level_value'] = dict(value=0, age_before=10, age_after=upper_bound)

    # # # model.vars += dismod3.ism.consistent(model)

    # # # dismod3.fit.fit_consistent(model, iter=iter, thin=thin, burn=burn)
    # # # bounds[upper_bound] = model

# # # # <codecell>

# # # figure(**book_graphics.half_page_params)

# # # param_list = [dict(age=100, linestyle='-', label='onset < 100'),
              # # # dict(age=65, linestyle='--', label='onset < 65'),
              # # # dict(age=45, linestyle=':', label='onset < 45')]
# # # for params in param_list:
    # # # model = bounds[params['age']]
    # # # p = dismod3.covariates.predict_for(model, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                       # # # 1., model.vars['p'], 0., 1.).T
    # # # i = dismod3.covariates.predict_for(model, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                       # # # 1., model.vars['i'], 0., 1.).T

    # # # subplot(1,2,1)
    # # # plot(arange(101), p.mean(axis=1), 'k', linestyle=params['linestyle'], linewidth=2, label=params['label'])
    
    # # # subplot(1,2,2)
    # # # plot(arange(101), i.mean(axis=1), 'k', linestyle=params['linestyle'], linewidth=2, label=params['label'])

# # # subplot(1,2,1)
# # # dismod3.graphics.plot_data_bars(model.get_data('p'))
# # # xlabel('Age (years)')
# # # ylabel('Prevalence (%)')
# # # yticks([0, .01, .02], [0, 1, 2])
# # # axis([-5,105,-.0005, .022])


# # # subplot(1,2,2)
# # # xlabel('Age (years)')
# # # ylabel('Incidence (per 10,000 PY)')
# # # yticks([0, .0005, .001, .0015], [0, 5, 10, 15])
# # # axis([-5,105,-.00001, .002])
# # # legend(loc='upper center', fancybox=True, shadow=True, bbox_to_anchor=(.8,1))

# # # # <markdowncell>

# # # # What happens if incidence has less smoothing?
# # # # --------------------------

# # # # <codecell>

# # # model = dismod3.data.load('/snfs1/Project/dismod/notebooks/models/bipolar')
# # # #model.keep(areas=['europe_western'], sexes=['male', 'total'], end_year=1997)


# # # model.parameters['i']['smoothness']['amount'] = 'Slightly'

# # # model.vars += dismod3.ism.consistent(model)

# # # dismod3.fit.fit_consistent(model, iter=iter, thin=thin, burn=burn)
# # # slight_smooth = model

# # # # <codecell>

# # # slight_p = dismod3.covariates.predict_for(slight_smooth, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                          # # # 1., model.vars['p'], 0., 1.).T
# # # mod_p = dismod3.covariates.predict_for(best, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                       # # # 1., best.vars['p'], 0., 1.).T

# # # plot(arange(101), slight_p.mean(axis=1), 'k-', linewidth=2, label='Slight Smoothing')
# # # plot(arange(101), mod_p.mean(axis=1), 'k--', linewidth=2, label='Moderate Smoothing')

# # # #hpd = mc.utils.hpd(slight_p - mod_p, .05)
# # # #plot(arange(101), hpd[:,0], 'k-', linewidth=1, label='95% HPD interval')
# # # #plot(arange(101), hpd[:,1], 'k-', linewidth=1)

# # # xlabel('Age (years)')
# # # ylabel('Prevalence (%)')
# # # yticks([0, .005, .01], [0, .5, 1])
# # # legend(loc='upper right', fancybox=True, shadow=True)
# # # #axis(params.get('axis', [-5,105,-.005,.06]))

# # # # <markdowncell>

# # # # Prevalence type
# # # # ---------------
# # # # 
# # # # To maximise data inclusion, both point and 12-month prevalence estimates were included. 
# # # # A study level covariate was used to create a crosswalk between point prevalence and 12-month
# # # # prevalence with the latter set as the desirable. The rationale for doing this relates to the
# # # # presentation of bipolar disorder as an episodic disease.   Cases fluctuate between manic and depressive episodes, 
# # # # interspersed by periods of residual symptoms {American Psychiatric Association, 2000 #31}. Although residual symptoms
# # # # can be less severe than manic and depressive episodes, they still lead to some disability and therefore burden.
# # # # However since an average episode of bipolar disorder is believed to last for about 3 months or more
# # # # {Angst, 2000 #196;Goodwin, 2008 #208}, estimates of point prevalence assessing symptoms within the past
# # # # month or less, will likely miss cases of bipolar disorder in a residual episode and underestimate prevalence. 
# # # # 
# # # # The inclusion of the prevalence type study level covariate had a noteworthy impact on the prevalence prior estimates, 
# # # # adjusting point prevalence estimates up by 20% (5%-33%) to ‘crosswalk’ to corresponding values of 12-month prevalence,
# # # # the specified gold standard (Figures 3 and 4). (Note-We now no longer see the effect of this covariate in the latest model)

# # # # <codecell>

# # # # should cv_past_year = 1 be the "reference" value? no

# # # #model.output_template['x_cv_past_year']=1

# # # # <codecell>

# # # 1.-exp(model.parameters['p']['fixed_effects']['x_cv_past_year']['mu'])

# <codecell>

model = dismod3.data.load('/home/j/Project/dismod/notebooks/models/bipolar')
df = model.get_data('p')

pl.figure(**book_graphics.half_page_params)
for i in range(2):
    pl.subplot(1,2,1+i)
    dismod3.graphics.plot_data_bars(df[df['x_cv_past_year'] == i])
    pl.xlabel('Age (years)', fontsize='xx-large')
    pl.ylabel('Prevalence (%)', fontsize='xx-large')
    pl.yticks([0, .01, .02, .03, .04], [0, 1, 2, 3, 4], fontsize='x-large')
    pl.axis([-5,105,-.0045, .045])
    if i == 0: subtitle('(a) Past-year prevalence')
    else: subtitle('(b) Past-month prevalence')
    #title('cv_past_year = %d'%i)

pl.subplots_adjust(wspace=.35, bottom=.14)    
pl.savefig('book/graphics/bipolar-data-by-cv.pdf')

# <codecell>

# # # model = dismod3.data.load('/snfs1/Project/dismod/notebooks/models/bipolar')

# # # # remove expert prior on pyp effect
# # # model.parameters['p']['fixed_effects'].pop('x_cv_past_year')

# # # model.vars += dismod3.ism.consistent(model)

# # # dismod3.fit.fit_consistent(model, iter=iter, thin=thin, burn=burn)

# # # # <codecell>

# # # no_py_ep = model

# <markdowncell>

# Compare alternative reference values
# -------------------------------------

# <codecell>

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

pl.figure(**book_graphics.half_page_params)
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

pl.subplot(1,2,2)
pl.xlabel('Age (years)', fontsize='xx-large')
pl.ylabel('Incidence \n (per 10,000 PY)'+'\n\n', ha='center',fontsize='xx-large')
pl.yticks([0, .0005, .001, .0015, .002], [0, 5, 10, 15, 20], fontsize='x-large')
pl.plot([-10],[-10],'k-', label='Past year')
pl.plot([-10],[-10],'k--', label='Past month')
pl.legend(bbox_to_anchor=(.42, 0, .5, .94), bbox_transform=pl.gcf().transFigure, shadow=True)
pl.axis([-5,105,-.00023, .0023])
pl.subplots_adjust(wspace=.35, bottom=.14)

pl.savefig('book/graphics/bipolar-ref-alts.pdf')

pl.show()

# # # <codecell>

# # model.vars['p']['beta'][1].__name__, model.vars['p']['beta'][1].stats()

# # # <codecell>

# # model = dismod3.data.load('models/bipolar')

# # model.vars += dismod3.ism.consistent(model)

# # dismod3.fit.fit_consistent(model, iter=iter, thin=thin, burn=burn)

# # # <codecell>

# # model.vars['p']['beta'][1].__name__, model.vars['p']['beta'][1].stats()

# # # <codecell>

# # py_ep = model

# # # <codecell>

# # model.parameters['p']['fixed_effects']

# # # <codecell>

# # no_py_ep.output_template['x_cv_past_year'] = 0.
# # py_ep.output_template['x_cv_past_year'] = 0.

# # # <codecell>

# # wo = dismod3.covariates.predict_for(no_py_ep, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                    # # 1., no_py_ep.vars['p'], 0., 1.).T
# # w = dismod3.covariates.predict_for(py_ep, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                   # # 1., py_ep.vars['p'], 0., 1.).T

# # plot(arange(101), wo.mean(axis=1), 'k-', linewidth=2, label='No prior on beta_py')

# # plot(arange(101), w.mean(axis=1), 'k--', linewidth=2, label='beta_py ~ Normal(-.45,.2)')

# # xlabel('Age (years)')
# # ylabel('Prevalence (%)')
# # yticks([0, .005, .01], [0, .5, 1])
# # legend(loc='upper right', fancybox=True, shadow=True)
# # #axis(params.get('axis', [-5,105,-.005,.06]))

# # # <codecell>

# # print 'With expert priors, pmp:pyp ratio is:', pl.exp(py_ep.vars['p']['beta'][1].stats()['mean']), pl.exp(py_ep.vars['p']['beta'][1].stats()['95% HPD interval'])

# # # <codecell>

# # print 'Without expert priors, pmp:pyp ratio is:', pl.exp(no_py_ep.vars['p']['beta'][1].stats()['mean']), pl.exp(no_py_ep.vars['p']['beta'][1].stats()['95% HPD interval'])

# # # <markdowncell>

# # # Lack of incidence data
# # # --------------
# # # 
# # # We found only 2 studies meeting our inclusion criteria for the annual
# # # incidence of bipolar disorder. Both these studies used samples from the United States
# # # of America and reported incidence estimates of 0% and 0.13% respectively {Lewinsohn, 1993 #321;Lewinsohn, 1995 #322}.
# # # Since we had more prevalence data in comparison, these 2 incidence studies were excluded from the modelling
# # # process and incidence was derived from prevalence, remission and mortality data. 
# # # 
# # # Figures 5 and 6 compare models where existing incidence estimates were used (model 17164) and where
# # # incidence was modelled using data from other parameters (model  16475). The difference between
# # # the empirical prior and posterior for incidence in figure 5 suggests a lack of consistency between the
# # # available incidence estimates and data from other parameters.  It also illustrates how DisMod can weight
# # # data from other parameters more heavily (i.e. prevalence) when deriving output for a parameter with
# # # very little data. The inclusion of existing incidence data also produced lower prevalence
# # # posteriors compared to when incidence was derived from other parameters (figure 6). 
# # # Theo/Jed I wasn’t too sure about this paragraph so some input here will be helpful. 

# # # <codecell>

# # model = dismod3.data.load('models/bipolar_orig')

# # # <codecell>

# # model.parameters = best.parameters

# # # <codecell>

# # model.vars += dismod3.ism.consistent(model)

# # %time dismod3.fit.fit_consistent(model, iter=iter, thin=thin, burn=burn)
# # w_incidence = model

# # # <codecell>

# # w_incidence.plot_asr('f')

# # # <codecell>

# # best.plot_asr('f')

# # # <codecell>

# # wo = dismod3.covariates.predict_for(best, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                    # # 1., best.vars['p'], 0., 1.).T
# # w = dismod3.covariates.predict_for(w_incidence, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                   # # 1., w_incidence.vars['p'], 0., 1.).T

# # plot(arange(101), wo.mean(axis=1), 'k-', linewidth=2, label='Without incidence data')

# # plot(arange(101), w.mean(axis=1), 'k--', linewidth=2, label='With incidence data')

# # xlabel('Age (years)')
# # ylabel('Prevalence (%)')
# # yticks([0, .005, .01], [0, .5, 1])
# # legend(loc='upper right', fancybox=True, shadow=True)
# # #axis(params.get('axis', [-5,105,-.005,.06]))

# # # <markdowncell>

# # # Lack of remission data
# # # -----------------------
# # # 
# # # The terms ‘residual’ and ‘remission’ have very different implications for GBD. As previously explained, 
# # # a residual state involves mild symptoms with mild disability which still contribute to burden. Remission is 
# # # equivalent to cure rather than a temporary reduction in symptom levels hence does not contribute to burden. 
# # # Since there is no consistent operationalisation of these terms in the bipolar literature, we were unable to include 
# # # any remission data in the bipolar modelling. Instead, expert guidance was sought to set a fix remission rate. Given 
# # # the chronicity of bipolar symptoms the model was trialed with remission set to 0 and to a maximum of 0.05 with the latter 
# # # yielding a more plausible fit to the existing data. Figure 7 illustrates the remission output with a maximum rate of 0.05 imposed.

# # # <codecell>

# # remission = {}

# # for upper_bound in [0, .05, .1]:
    # # model = dismod3.data.load('models/bipolar')
    # # #model.keep(areas=['europe_western'], sexes=['male', 'total'], end_year=1997)

    # # model.parameters['r']['level_bounds'] = dict(lower=0, upper=upper_bound)

    # # model.vars += dismod3.ism.consistent(model)

    # # dismod3.fit.fit_consistent(model, iter=iter, thin=thin, burn=burn)
    # # remission[upper_bound] = model

# # # <codecell>

# # figure(**book_graphics.half_page_params)

# # param_list = [dict(ub=.1, linestyle='-', label='remission < 10 '),
              # # dict(ub=.05, linestyle='--', label='remission < 5'),
              # # dict(ub=0, linestyle=':', label='remission = 0')]
# # for params in param_list:
    # # model = remission[params['ub']]
    # # p = dismod3.covariates.predict_for(model, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                       # # 1., model.vars['p'], 0., 1.).T
    # # r = dismod3.covariates.predict_for(model, 'all', 'total', 'all', 'europe_western', 'male', 1990,
                                       # # 1., model.vars['r'], 0., 1.).T

    # # subplot(1,2,1)
    # # plot(arange(101), p.mean(axis=1), 'k', linestyle=params['linestyle'], linewidth=2, label=params['label'])
    
    # # subplot(1,2,2)
    # # plot(arange(101), r.mean(axis=1), 'k', linestyle=params['linestyle'], linewidth=2, label=params['label'])

# # subplot(1,2,1)
# # dismod3.graphics.plot_data_bars(model.get_data('p'))
# # xlabel('Age (years)', fontsize='xx-large')
# # ylabel('Prevalence (%)', fontsize='xx-large')
# # yticks([0, .005, .01, .015], [0, .5, 1, 1.5], fontsize='x-large')
# # axis([-5,105,-.0005, .02])
# # legend(loc='upper center', fancybox=True, shadow=True, bbox_to_anchor=(.8,1))

# # subplot(1,2,2)
# # xlabel('Age (years)', fontsize='xx-large')
# # ylabel('Remission (per 100 PY)', fontsize='xx-large')
# # yticks([0, .05, .1], [0, 5, 10], fontsize='xx-large')
# # #axis([-5,105,-.001, .12])

# # # <markdowncell>

# # # Global application of excess mortality data
# # # -------------------------------------------
# # # 
# # # The excess mortality data available for bipolar disorder were standard mortality ratios (SMRs) 
# # # from developed parts of the world.  SMRs are derived from mortality risks in cohort studies of 
# # # patients with bipolar disorder expressed as a ratio to expected deaths if general population 
# # # mortality risks would have applied.  Due to lack of data, the available SMRs were applied equally 
# # # across all GBD regions.  Although this facilitated the modelling process in the absence of more 
# # # representative data, it led to excessively high rates of excess mortality in Sub-Saharan Africa 
# # # due to high mortality rates from HIV/AIDS in young and middle aged adults.
# # # 
# # # To deal with this, we first ran a model incorporating the available SMRs. The region-, sex-, year- 
# # # and age-specific excess mortality estimates DisMod generated based on this data were then exported 
# # # and re-inserted into the model as data points to replace the SMR data.  At this stage, estimates 
# # # from Sub-Saharan Africa Central, Southern and East were replaced with estimates from Sub-Saharan Africa West 
# # # where the effect brought about by mortality due to HIV/AIDS was not as strong. Figures 8 and 9 show prevalence 
# # # posteriors with the original SMR data used and with the derived excess mortality estimates used respectively.  
# # # Theo is this better?

# # # <codecell>

# # # TODO: plot comparing prevalence in SSA w these two approaches

# # # <markdowncell>

# # # Time trend on prevalence
# # # ------------------------
# # # 
# # # After further interrogation of the model, an additional assumption was made.  
# # # An expert prior specifying no time difference was included as experts were 
# # # not confident that differences observed over time were a true time trend or 
# # # due to differences in study methodology over time. Figures 10 to 13 show prevalence 
# # # posteriors with and without the time effect covariate. Although this effect seems 
# # # small, it potentially has significant consequences at a population level and on the 
# # # final estimates of disease burden for 1990 and 2005.
# # # 
# # # An alternative approach here would have been to control for methodological differences 
# # # between studies through the inclusion of other study-level covariates which adjust for 
# # # study quality. For other mental disorders the choice of study quality covariates was based 
# # # on analysis of variance outside of dismod, highlighting variables (e.g. sample 
# # # representativeness, diagnostic criteria) with a significant impact on prevalence.  
# # # Given the lack of available data for bipolar disorder, analyses outside of DisMod provided v
# # # ery little rationale for the selection of covariates {Ferrari, 2010 #404}.   

# # # <codecell>

# # model = best
# # df = model.vars['p']['data']
# # df['mu_pred'] = model.vars['p']['p_pred'].stats()['mean']
# # df['err'] = df['value'] - df['mu_pred']

# # # <codecell>

# # df['value'] = model.input_data['value']

# # # <codecell>

# # plot(df['year_start'].__array__(), df['err'].__array__(), 'ks', mec='w')
# # hlines([0], 1975, 2010, linestyle='--')
# # axis([1975,2010,-.016,.025])
# # xlabel('Start Year of Observation')
# # ylabel('Residual (observed - predicted)')

# # # <markdowncell>

# # # Outlier covariate
# # # -----------------
# # # 
# # # What is the story here?

# # # <codecell>

# # model = dismod3.data.load('models/bipolar')

# # # <codecell>

# # dismod3.graphics.scatter(model.get_data('p'), 'age_start', 'value', groupby='x_cv_outlier')

# # # <codecell>

# # # TODO: compare model with and without outlier covariate effect, with generalized NB instead of shift

# # # <markdowncell>

# # # Additional covariates to investigate?

# # # <markdowncell>

# # # Global heterogeneity
# # # --------------------

# # # <codecell>

# # heterogeneity = {'Moderately': best}

# # for het in ['Slightly', 'Very']:
    # # model = dismod3.data.load('models/bipolar')
    # # model.parameters['p']['heterogeneity'] = het

    # # model.vars += dismod3.ism.consistent(model)

    # # dismod3.fit.fit_consistent(model, iter=iter, thin=thin, burn=burn)
    # # heterogeneity[het] = model

# # # <codecell>

# # figure(figsize=(17,8.5))
# # for i, het in enumerate(['Slightly', 'Moderately', 'Very']):
    # # subplot(1,3,i+1)
    # # model = heterogeneity[het]
    # # for sr in model.hierarchy['all']:
        # # pred = dismod3.covariates.predict_for(model, 'all', 'total', 'all', sr, 'male', 1990,
                                              # # 1., model.vars['p'], 0., 1.).T

        # # plot(arange(101), pred.mean(axis=1), '-', linewidth=2, label=sr)

    # # title(het)
    # # xlabel('Age (years)')
    # # ylabel('Prevalence (%)')
    # # yticks([0, .005, .01], [0, .5, 1])
    
    # # axis(params.get('axis', [-5,105,-.005,.025]))
    
# # legend(loc='upper right', fancybox=True, shadow=True)

# # # <codecell>


