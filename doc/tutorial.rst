Tutorial with motivating example
================================

The goal of this document it give a concise demonstration of the 
strengths and limitations of DisMod-MR, the descriptive
epidemiological meta-regression tool developed for the Global Burden of Disease,
Injuries, and Risk Factors 2010 (GBD 2010) Study.

Descriptive epidemiological meta-regression of Parkinson's Disease
------------------------------------------------------------------

A systematic review of Parkinson's disease (PD) was conducted as part of the GBD 2010
Study. The results of this
review---data on the prevalence, incidence, and standardized mortality ratio of
PD---needed to be combined to produce estimates of disease prevalence by
region, age, sex, and year.  These prevalence estimates were combined
with disability weights to measure years lived with disability (YLDs),
which were then combined with estimates of years of life lost (YLLs)
to produce estimates of the burden of PD quantified in disability-adjusted life-years (DALYs).

PD is a neurodegenerative disorder that includes symptoms of motor
dysfunction, such as tremors, rigidity, and akinesia, in the early
stages of the disease.  As the disease develops, most patients also
develop nonmotor symptoms, such as cognitive decline, dementia,
autonomic failure, and disordered sleep-wake regulation.  The standard
definition for PD diagnosis includes at least two of four cardinal
signs---resting tremor, bradykinesia, rigidity, and postural abnormalities.
There is no cure or treatments to slow the progression of the disease;
however, motor symptoms and disability may be improved with
symptomatic therapy.

.. sourcecode:: python

	In [1]: # setup dismod3 environment
		# TODO: make this automatic

		%cd /homes/abie/gbd_dev/gbd
		import sys
		sys.path += ['book']

		import pymc as mc, dismod3, book_graphics, pandas
		
		mpl.rcParams['axes.titlesize'] = 'xx-large'; mpl.rcParams['axes.labelsize'] = 'xx-large'; 
		mpl.rcParams['xtick.labelsize'] = 'x-large'; mpl.rcParams['ytick.labelsize'] = 'x-large'; 
		mpl.rcParams['legend.fancybox'] = True; mpl.rcParams['legend.fontsize'] = 'large'; mpl.rcParams['font.size'] = 12
		
		def my_axis(ymax): axis([-5,105,-ymax/10.,ymax])

		def subtitle(s):
			""" title where the panel names appear within each panel"""
			l,r,b,t=axis()
			x = l + (r-l)*.05
			y = t - (t-b)*.05
			text(x, y, s, ha='left', va='top', size='x-large')

DisMod-MR uses the integrative systems modeling (ISM) approach to produce simultaneous
estimates of disease incidence, prevalence, remission, and mortality. The hallmark of
ISM is incorporating all available data.  In the case of Parkinson's Disease this
consists of population level measurements of incidence, prevalence, standardized mortality rate (SMR),
and cause-specific mortality rate (CSMR).

I will begin with a look at a subset of this data, however.  Only that from females in the Europe, Western GBD region.

.. sourcecode:: python

	In [2]: model = dismod3.data.load('/home/j/Project/dismod/output/dm-40552')
		model.keep(areas=['europe_western'], sexes=['female'])

Of the 235 rows of data, here is how the values breakdown by data type:

.. sourcecode:: python

	In [3]: summary = model.input_data.groupby('data_type')['value'].describe()
		round_(summary,3).sort('count', ascending=False)
	Out [3]:       count  mean   std    min   10%    50%    90%    max  
		m_all  1848   0.076  0.134  0     0.001  0.01   0.279  0.554
		csmr   1638   0      0      0     0      0      0      0.002
		p      660    0.01   0.013  0     0      0.004  0.028  0.071
		i      99     0.001  0.001  0     0      0      0.002  0.005
		smr    13     2.842  1.188  1.06  1.32   2.55   3.88   5.41

More than half of the available data for this region is prevalence data.  I'll take a closer look at that now.

.. sourcecode:: python

	In [4]: groups = model.get_data('p').groupby('area')
		round_(groups['value'].describe(),3).sort('50%', ascending=False)
	Out [4]:                count  mean   std    min    10%    50%    90%    max  
		NLD             12     0.03   0.02   0.002  0.006  0.034  0.05   0.053
		europe_western  6      0.02   0.01   0.006  0.008  0.024  0.028  0.029
		FRA             7      0.021  0.024  0.002  0.002  0.018  0.042  0.071
		ESP             29     0.018  0.016  0      0      0.013  0.04   0.05 
		ITA             39     0.007  0.01   0      0      0.003  0.022  0.041
		DEU             1      0.003  NaN    0.003  0.003  0.003  0.003  0.003
		GBR             32     0.004  0.006  0      0      0.002  0.011  0.021
		FIN             1      0.002  NaN    0.002  0.002  0.002  0.002  0.002
		NOR             1      0.001  NaN    0.001  0.001  0.001  0.001  0.001
		PRT             11     0.002  0.003  0      0      0      0.006  0.009

There is a wide range in median values, which reflects a combination of country-to-country variation
and compositional bias.  I'll compare data from NLD and PRT, countries with a moderate amount of data each,
that have median measured values of 33 per 1000 and 5 per 100,000 respectively.

.. sourcecode:: python

	In [5]: NLD = groups.get_group('NLD')
		PRT = groups.get_group('PRT')
		
	In [6]: figure(**book_graphics.half_page_params)

		for i, s, d in [[1, '(a) NLD', NLD], [2, '(b) PRT', PRT]]:
			subplot(1,2,i)
			dismod3.graphics.plot_data_bars(d)
			xlabel('Age (years)')
			ylabel('Prevalence (%)')
			yticks([0, .02, .04, .06, .08], [0, 2, 4, 6, 8])
			my_axis(.06)
			subtitle(s)
			grid()
	Out [6]: 

.. image:: data_nld_prt.pdf
	:align: center

In these plots, every row of data collected in systematic review is represented as a two squares
joined by a horizontal line.  The distance above the $x$-axis denotes the measured prevalence, and the
left and right point shows the starting and ending age of the age group observed.

These plots show why the medians of the measurements are so different: the age groups reported for PRT 
include observations of the very young, where PD is not present.  The age groups reported for NLD focus
particularly on older ages.  For the age groups that match, the levels are quite similar.

A model for age-specific parameters when measurements have heterogeneous age groups
-----------------------------------------------------------------------------------

DisMod-MR has four features that make it particularly suited for estimating age-specific prevalence of PD from this data:

- Piecewise linear spline model for change in prevalence as a function of age
- Age-standardizing model of age-group heterogeneity represents the heterogeneous age groups collected in systematic review
- Country-level random effects for true variation in prevalence between countries
- Negative binomial model of data, which provides data-driven estimation of non-sampling error in measurements
  and elegantly handles measurements of 0

I will now fit the prevalence data with DisMod-MR's age-standardizing negative binomial random effect spline model and
compare the estimates to the observed data.  Then I will use the results of the fit model to explore the four features listed above.

.. sourcecode:: python

	In [7]: # remove fixed effects for this example, I will return to them below
		model.input_data = model.input_data.filter(regex='(?!x_)')
	
	In [8]: model.vars += dismod3.ism.age_specific_rate(model, 'p')
		%time dismod3.fit.fit_asr(model, 'p')
	
	In [9]: # plot age-specific prevalence estimates over data bars
		figure(**book_graphics.half_page_params)

		dismod3.graphics.plot_data_bars(model.get_data('p'), color='grey')
		pred = dismod3.covariates.predict_for(model, model.parameters['p'], 'all', 'female', 2005,
							'europe_western', 'female', 2005, 1.,
							model.vars['p'], 0., 1.)    # TODO: simplify this method!
		hpd = mc.utils.hpd(pred, .05)

		plot(arange(101), pred.mean(axis=0), 'k-', linewidth=2, label='Posterior Mean')
		plot(arange(101), hpd[:,0], 'k--', linewidth=1, label='95% HPD interval')
		plot(arange(101), hpd[:,1], 'k--', linewidth=1)

		xlabel('Age (years)')
		ylabel('Prevalence (%)\n', ha='center')
		yticks([0, .02, .04, .06, .08], [0, 2, 4, 6, 8])
		my_axis(.06)
		grid()
		legend(loc='upper left')
	Out [9]: TODO: PLOT HERE!
	
	In [10]: p_only = model  # store results for future comparison

This estimate shows the nonlinear increase in prevalence as a function of age, where the slope of the
curve increases at age 60.  A nonlinear estimate like this is possible thanks to DisMod-MR's piecewise linear
spline model.

The age-standardizing model for heterogeneous age groups is also important for
such settings; a naive approach, such as using the age interval midpoint, would result in under-estimating
the prevalence for age groups that include both individuals older and younger than 60.

The exact age where the slope of the curve changes is _not_ entirely data driven in this example.  The knots
in the piecewise linear spline model were chosen a priori, on the following grid:	
	
.. sourcecode:: python

	In [11]:  model.parameters['p']['parameter_age_mesh']
	Out [11]: [0, 30, 45, 60, 80, 100]

A sparse grid allows faster computation, but a dense grid allows more expressive age pattens.  Choosing
the proper balance is one challenge of a DisMod-MR analysis.  This is especially true for sparse,
noisy data, where too many knots allow the model to follow noisy idiosyncrasies of the data.  DisMod-MR
allows for penalized spline regression to help with this choice.

The country-level random effects in this model capture country-to-country variation in PD prevalence.
This variation is not visible in the graphic above, which shows the regional aggregation of country-level
estimates (using a population weighted average that takes uncertainty into account).

The country-level random effects take the form of intercept shifts in log-prevalence space, with values
showing in the following:

.. sourcecode:: python

	In [12]: df = pandas.DataFrame(index=[alpha_i.__name__ for alpha_i in model.vars['p']['alpha']],
                      columns=['mean', 'lb', 'ub'])
		 for alpha_i in model.vars['p']['alpha']:
		      stats = alpha_i.stats()
		      df.ix[alpha_i.__name__] = (stats['mean'], stats['95% HPD interval'][0], stats['95% HPD interval'][1])
	
	In [13]: round_(df,1).sort('mean', ascending=False)
	Out [13]:             mean  lb   ub 
		 alpha_p_NLD  0.1  -0.1  0.4
		 alpha_p_NOR  0    -0.2  0.2
		 alpha_p_ITA  0    -0.1  0.2
		 alpha_p_FIN  0    -0.2  0.2
		 alpha_p_FRA  0    -0.2  0.1
		 alpha_p_PRT  0    -0.2  0.3
		 alpha_p_ESP  0    -0.2  0.2
		 alpha_p_DEU  0    -0.2  0.2
		 alpha_p_GBR -0.1  -0.3  0.1

This shows that although none of the country-to-country variation is significant, NLD does appear 10% higher than the regional
average.  PRT appears no lower than average, however, and the difference in raw data medians was indeed due to compositional
bias.  The country with lowest estimate is GBR, estimated to be 10% lower than average.

The fourth feature of the model which I want to draw attention to here is the negative binomial model of data,
which deals with measurements of zero prevalence in a principled way.  Prevalence studies are reporting transformations
of count data, and count data can be zero.  In the case of prevalence of PD in 30- to 40-year-olds, it often _will_ be zero.

.. sourcecode:: python

	In [14]: model.get_data('p').sort('age_start').filter(['age_start', 'age_end', 'area', 'value']).head(15)
	Out [14]:      age_start  age_end  area  value   
		  276  0          49       ITA   0       
		  371  0          49       ITA   0       
		  394  0          54       ITA   0       
		  559  0          4        PRT   0       
		  563  5          9        PRT   0       
		  574  10         14       PRT   0       
		  575  15         24       PRT   0       
		  558  25         34       PRT   5e-05   
		  68   30         39       ESP   7.1e-05 
		  129  30         99       FIN   0.00189 
		  183  30         39       GBR   0       
		  207  30         50       GBR   3e-05   
		  245  30         40       GBR   1e-05   
		  311  30         99       ITA   0.002203
		  313  30         99       ITA   0.002183

The negative binomial model has an appropriately skewed distribution, where prevalence measurements 
of zero are possible, but measurements of less than zero are not possible.  To demonstrate how this
functions, the next figure shows the "posterior predictive distribution" for the measurements above,
i.e. sample values that the model predicts would be found of the studies were conducted again under
the same conditions.

.. sourcecode:: python

	In [15]: pred = model.vars['p']['p_pred'].trace()
		 obs = array(model.vars['p']['p_obs'].value)
		 ess = array(model.vars['p']['p_obs'].parents['n'])

	In [16]: figure(**book_graphics.half_page_params)

		 sorted_indices = obs.argsort().argsort()
		 jitter = mc.rnormal(0, .1**-2, len(pred))

		 for i,s_i in enumerate(sorted_indices):
			 plot(s_i+jitter, pred[:, i], 'ko', mew=0, alpha=.25, zorder=-99)

		 errorbar(sorted_indices, obs, yerr=1.96*sqrt(obs*(1-obs)/ess), fmt='ks', mew=1, mec='white', ms=5)

		 xticks([])
		 xlabel('Measurement')
		 ylabel('Prevalence (%)\n', ha='center')
		 yticks([0, .02, .04, .06, .08], [0, 2, 4, 6, 8])
		 axis([25.5,55.5,-.01,.1])
		 grid()
		 subtitle('Posterior Predictive distribution')
	Out [16]: TODO: PLOT HERE!
	
Additional features of DisMod-MR
--------------------------------

Four additional features of DisMod-MR that are important for many settings are:

- informative priors
- fixed effects to cross-walk between different studies
- fixed effects to predict out of sample
- fixed effects to explain the level of variation

Informative priors are useful for modeling disease with less data available than PD, for example to include
information that prevalence is zero for youngest ages, or than prevalence must be increasing as a function of
age between certain ages.

The informative priors are also key to the "empirical Bayes" approach to modeling age-specific differences between
difference GBD regions.  In this setting, a model using all the world's data is used to produce estimates for each region,
and these estimates are used as priors in region-specific models together with the data relevant to that region only.

"Cross-walk" fixed effects can correct for biases introduced by multiple outcome measures.  For example, in the PD dataset,

.. sourcecode:: python

	In [17]: model = dismod3.data.load('/home/j/Project/dismod/output/dm-40552')
	
	In [18]: crosswalks = list(model.input_data.filter(like='x_cv').columns)
		 groups = model.get_data('p').groupby(crosswalks)

	In [19]: crosswalks
	
	In [20]: round_(groups['value'].describe(),3)
	Out [20]:          count  mean   std    min    10%    50%    90%    max  
		  0  0  0  435    0.011  0.013  0      0      0.005  0.03   0.071
		  0  1  0  22     0.006  0.004  0      0      0.007  0.01   0.012
		  1  0  0  64     0.01   0.016  0      0      0.004  0.021  0.071
		  1  1  0  138    0.007  0.011  0      0      0.002  0.022  0.06 
			    1  1      0.005  NaN    0.005  0.005  0.005  0.005  0.005

Predictive fixed effects attempt to relate levels of disease to known covariates, for example
caffeine consumption per capita or smoking prevalence.

.. sourcecode:: python

	In [21]: figure(**book_graphics.half_page_params)
		 data = model.get_data('p')

		 for i, s, cv in [[1, '(a) FAO Stimulants', 'x_ihme_fao_stimulants_kcal_26oct11'], [2, '(b) Smoking Prev', 'x_smoking_prev']]:
		 	subplot(1,2,i)
		 	plot(array(data[cv]), array(data['value']), 'ks')
			
		 	xlabel(s)
		 	ylabel('Prevalence (%)')
		 	yticks([0, .02, .04, .06, .08], [0, 2, 4, 6, 8])
		 	xticks([])
		 	subtitle(s)
		 	grid()
	Out [21]: TODO: PLOT HERE!
	
As the scatter shows, the relationships are not very strong in this case.

Variation fixed effects function similarly, but are used to predict the negative binomial overdispersion
instead of the bias in the negative binomial mean.

Incorporating data on parameters other than prevalence
------------------------------------------------------

So far this example has focused on modeling the prevalence of PD from the
prevalence data alone.  However, this represents about half of the available
data.  There is also information on incidence, SMR, and CSMR, which has not
yet been incorporated.

DisMod-MR is capable of including all of the available data, using a compartmental
model of disease moving through a population.  This model formalizes the observation
that prevalent cases must once have been incident cases, and continue to be prevalent
cases until remission or death.

In this model, incidence, remission, and excess-mortality are age-standardizing negative binomial random effect spline models,
while prevalence, SMR, CSMR, and other parameters come from the solution to a system of ordinary differential equations.

The results of this model are smoother prevalence curves that take longer to calculate.

.. sourcecode:: python

	In [22]: figure(**book_graphics.full_page_params)
		 subplot(2,2,1); dismod3.graphics.plot_data_bars(model.get_data('p')); xlabel('Age (years)'); ylabel('Prevalence (%)'); yticks([0, .02, .04, .06, .08], [0, 2, 4, 6, 8]); my_axis(.09); subtitle('(a)'); grid()
		 subplot(2,2,2); dismod3.graphics.plot_data_bars(model.get_data('i')); xlabel('Age (years)'); ylabel('Incidence \n(per 10,000 PY)\n\n', ha='center'); yticks([0, .0013, .0026, .0039, .0052], [0, 13, 26, 39, 52]); my_axis(.0055); subtitle('(b)'); grid()
		 subplot(2,2,3); dismod3.graphics.plot_data_bars(model.get_data('csmr')); xlabel('Age (years)'); ylabel('Cause-specific mortality \n(per 10,000 PY)\n\n', ha='center'); yticks([0, .0004, .0008, .0012, .0016], [0, 4, 8, 12, 16]); my_axis(.0018); subtitle('(c)'); grid()
		 subplot(2,2,4); dismod3.graphics.plot_data_bars(model.get_data('smr')); xlabel('Age (years)'); ylabel('Standardized \nmortality ratio\n\n', ha='center'); yticks([0, 1, 2, 3, 4]); my_axis(4.5); subtitle('(d)'); subplots_adjust(hspace=.35,wspace=.35); grid()
	Out [22]: TODO: PLOT HERE!
	
	In [23]: model.input_data = model.input_data.filter(regex='(?!x_)')
		 model.vars += dismod3.ism.consistent(model)
		 %time dismod3.fit.fit_consistent(model)

	In [24]: figure(**book_graphics.full_page_params)

		 param_list = [dict(type='p', title='(a)', ylabel='Prevalence (%)', yticks=([0, .01, .02], [0, 1, 2]), axis=[30,101,-0.001,.025]),
			  dict(type='i', title='(b)', ylabel='Incidence \n(per 1000 PY)', yticks=([0, .001,.002, .003, .004], [0, 1, 2, 3, 4]), axis=[30,104,-.0003,.0055]),
			  dict(type='pf', title='(c)', ylabel='Cause-specific mortality \n(per 1000 PY)', yticks=([0, .001,.002], [0, 1, 2]), axis=[30,104,-.0002,.003]),
			  dict(type='smr', title='(d)', ylabel='Standardized \nmortality ratio', yticks=([1, 2, 3,4, ], [1, 2,3, 4]), axis=[35,104,.3,4.5]),
			  ]

		 for i, params in enumerate(param_list):
		 	ax = subplot(2,2,i+1)
		 	if params['type'] == 'pf': dismod3.graphics.plot_data_bars(model.get_data('csmr'), color='grey')
		 	else: dismod3.graphics.plot_data_bars(model.get_data(params['type']), color='grey')
			
			 if params['type'] == 'smr': model.pred = dismod3.covariates.predict_for(model, model.parameters.get('smr', {}), 'all', 'female', 2005, 
												   'europe_western', 'female', 2005, 1.,  model.vars['smr'], 0., 100.).T
			 else : model.pred = dismod3.covariates.predict_for(model, model.parameters[params['type']], 'all', 'female', 2005, 
									    'europe_western', 'female', 2005, 1.,  model.vars[params['type']], 0., 1.).T

			 plot(arange(101), model.pred.mean(axis=1), 'k-', linewidth=2, label='Posterior Mean')
			 hpd = mc.utils.hpd(model.pred.T, .05)
			 plot(arange(101), hpd[:,0], 'k-', linewidth=1, label='95% HPD interval')
			 plot(arange(101), hpd[:,1], 'k-', linewidth=1)

			 xlabel('Age (years)')
			 ylabel(params['ylabel']+'\n\n', ha='center')
			 axis(params.get('axis', [-5,105,-.005,.06]))
			 yticks(*params.get('yticks', ([0, .025, .05], [0, 2.5, 5])))
			 subtitle(params['title'])
			 grid()
			
		 subplots_adjust(hspace=.35, wspace=.35)
	Out [24]: TODO: PLOT HERE!
	
	In [25]: p_with = model
	
The most notable difference between the estimates from this model and from the model
that used prevalence data only is that this model produces estimates of incidence and
mortality in addition to prevalence.  In many cases, the model also produces estimates
of the remission rate as well, but there is no remission of PD, so the estimates of zero
are not very interesting in this example.  It is another place that informative priors are useful,
however.

There are also differences between the means and uncertainty intervals estimated by these methods,
which show that the additional data is important.  Although the prevalence data alone predicts 
age-specific prevalence that peaks at 2%, when the incidence and mortality data is also included, the
maximum prevalence is a bit lower, closer to 1.5%.

.. sourcecode:: python

	In [26]: p1 = dismod3.covariates.predict_for(p_only, model.parameters['p'], 'all', 'female', 2005, 
							'europe_western', 'female', 2005, 1., p_only.vars['p'], 0., 1.)

		 p2 = dismod3.covariates.predict_for(p_with, model.parameters['p'], 'all', 'female', 2005, 
							'europe_western', 'female', 2005, 1., p_with.vars['p'], 0., 1.)

	In [27]: plot(p1.mean(axis=0), 'k--', linewidth=2, label='Only prevalence')
		 plot(p2.mean(axis=0), 'k-', linewidth=2, label='All available')

		 xlabel('Age (years)')
		 ylabel('Prevalence (%)\n\n', ha='center')
		 yticks([0, .01, .02], [0, 1, 2])
		 axis([30,101,-0.001,.025])
		 legend(loc='upper left')
		 grid()
	Out [27]: TODO: PLOT HERE!
	
Because the data is so noisy, the differences between the mean estimates of these different models are not significant; the posterior distributions
have considerable overlap.  At age 80, for example, the posterior distributions for age-80 prevalence are estimated as the following:

.. sourcecode:: python

	In [28]: hist(100*p1[:,80], normed=True, histtype='step', label='Only prevalence', linewidth=3, color=array([239., 138., 98., 256.])/256)
		 hist(100*p2[:,80], normed=True, histtype='step', label='All available', linewidth=3, color=array([103, 169, 207, 256.])/256)
		 title('PD prevalence at age 80')
		 xlabel('Prevalence (%)\n\n', ha='center')
		 ylabel('Probability Density')
		 legend(loc='upper right')
		 grid()
	Out [28]: TODO: PLOT HERE!
	
Conclusion
----------

I hope that this example is a quick way to see the strengths and weaknesses of DisMod-MR.
This model is particularly suited for estimating descriptive epidemiology of diseases
with sparse, noisy data from multiple, incompatible sources.

I am currently working to make it faster, as well as to improve the capabilities for modeling
changes between regions over time.