High-level functions
====================

data.py
-------

.. automodule:: data

.. autoclass:: ModelData
	:members:
	:exclude-members: from_gbd_json, from_gbd_jsons

.. _hierarchy-label:

hierarchy.json 
**************

THIS FILE CONFUSES ME!!!!

The file contains a dictionary with keys 'nodes' and 'edges'.

'nodes' is a dictionary of each item in hierarchy.  Each item in the hierarchy is a dictionary.  Countries contain the `pop`, the country population. Regions and superregions are empty dictionaries.
'edges' is a list of

**Example:**

.. sourcecode:: python

	{"nodes": [], "edges": []}

.. _input_data-label:

input_data.csv 
**************

:Required fields:
  - `Region` : enter 'all' if parameter is to apply to all regions
  - `Country ISO3 Code` : ISO3 country code, enter 'all' if parameter is to apply to all regions
  - `GBD Cause` : open text field, enter the name of the cause
  - `Year Start` : year when data collection started, [1950, 2010], `Year Start` MUST <= `Year End`
  - `Year End` : year when data collection ended, [1950, 2010], `Year End` MUST >= `Year Start`
  - `Parameter` : one of 'Prevalence', 'Incidence', 'Remission', 'Mortality', 'RR' (relative risk), 'SMR' (standardized mortality ratio), 'With Condition Mortality', 'Duration', 'cause-specific mortality', or 'Excess mortality'
  - `Age Start` : beginning of age group, [0, 100], `Age Start` MUST <= `Age End`
  - `Age End` : end of age group, [0, 100], `Age End` MUST >= `Age Start`
  - `Sex` : one of 'Female', 'Male', or 'Total' (no gender differentiation)
  - `Parameter Value` : value of parameter, preferable to express 'per 1', MUST be >=0
  - `Units` : units of parameter value, MUST be >=1
  - `Lower CI` : value of lower confidence interval, if not stated leave blank, MUST be > 0 <= `Parameter Value`
  - `Upper CI` : value of upper confidence interval, if not stated leave blank, MUST be >= `Parameter Value`
  - `Standard error` : value of standard error, preferable to express 'per 1'
  
:Optional fields:
  - `Study ID` : article id, must be >= 0
  - `Country` : name of country
  - `Year of publication` : year when article was published (yyyy)
  - `Citation` : study reference details such as authors, title, journal name, year, volume no., issue no., and page no.
  - `Parameter type` : example include: point, past month, past year, lifetime, etc.
  - `Ignore` : binary, enter 1 if row is to be ignored
  - `%CI` : level of confidence interval, e.g., 90, 95, etc.
  - `Effective Sample Size` : effective sample size
  - `Denominator` : age-sex-specific sample size for age-sex band, i.e., N
  - `Numerator` : number of cases for the specific age-sex band with the condition
  - `Standardized?` : either 'yes' or 'no', indicates if the reported resulted were standardized
  - `Design Factor` : design factor, if not reported leave blank
  - `Study type` : type of study, e.g., survey, etc.
  - `Data ascertainment` : location of sampling population, community sample, clinic, etc.
  - `Overall study sample size` : number of people from whom data were collected
  - `Case Definition` : description of case definition
  - `Study population description` : breif description of study population, i.e. age range, geographical location, etc.
  - `Coverage string` : area of coverage, e.g., national, sub-national, regional, community, etc.
  - `Urbanicity string` : urbanicity of `Coverage string`, e.g., urban, rural, mixed
  - `Person-years of cases` : value of person-years of the condition (exposed population), additional column for mortality
  - `Person-years of population` : value of person-years of the population(non-exposed population), additional column for mortality
  - `Observed deaths` : number of deaths, additional column for mortality
  - `Follow-up period` : value of follow-up perion in years, additional column for mortality or remission
  - `Remitted %` : percentage remitted, additional column for remission
  - `Remission per year` : proportion calculated as -LN(1-(`Remitted %`/100))/`Follow-up period`, additional column for remission

.. note::
  - The `GBD Cause` name MUST be exactly the same for all rows.
  - If there is no upper limit for age, please enter 99, e.g., if age groups are 15-49 and 50+, use 15-49 and 50-99.
  - There MUST be an estimate of uncertainty.  If standard error and confidence intervals are not reported, enter an effective sample size.  If no effective sample size is reported, use the lowest effective sample size reported in other studies as a conservative estimate.

.. _nodes_to_fit-label:

nodes_to_fit.json
*****************

The file contains a list of areas contained in the hierarchy at which results are returned.
  
**Example:**

.. sourcecode:: python

	["all",  "north_america_high_income", "australasia", "europe_western", "asia_pacific_high_income", 
	 "latin_america_southern", "europe_central", "europe_eastern", "asia_central", "sub-saharan_africa_central", 
	 "sub-saharan_africa_east", "sub-saharan_africa_southern", "sub-saharan_africa_west", 
	 "north_africa_middle_east", "asia_south", "asia_southeast", "asia_east", "oceania", 
	 "latin_america_andean", "latin_america_central", "latin_america_tropical", "caribbean"]

.. _output_template-label:

output_template.csv 
*******************

This file contains the same headings as :ref:`input_data-label` but the contents are the appropriate reference values.

.. _parameters-label:

parameters.json 
***************

The file contains a dictionary with keys corresponding to disease parameters (such as 'rr', 'f', 'i', 'p', 'r', 'pf', and 'X').  Each key is a dictionary which contains priors for that parameter.  There is also a key 'ages', which contains a list of all ages. 

**Example:**

.. sourcecode:: python

	{
	  "p": {
		"increasing": {"age_start": 0, "age_end": 0}, 
		"fixed_effects": {}, 
		"level_bounds": {"upper": 1.0, "lower": 0.0}, 
		"y_maximum": 1.0, 
		"level_value": {"age_after": 100, "value": "0.0", "age_before": 1}, 
		"random_effects": {}, 
		"decreasing": {"age_start": 0, "age_end": 0}, 
		"parameter_age_mesh": [0, 1, 5, 10, 15, 20, 25, 35, 45, 55, 65, 75, 85, 100], 
		"heterogeneity": "Slightly", 
		"smoothness": {"age_start": 0, "amount": "Slightly", "age_end": 100}
	  },
	  "i": {
		"increasing": {"age_start": 0, "age_end": 0}, 
		"fixed_effects": {}, 
		"level_bounds": {"upper": 1.0, "lower": 0.0}, 
		"y_maximum": 1.0, 
		"level_value": {"age_after": 100, "value": "0.0", "age_before": 0}, 
		"random_effects": {}, 
		"decreasing": {"age_start": 0, "age_end": 0}, 
		"parameter_age_mesh": [0, 1, 5, 10, 15, 20, 25, 35, 45, 55, 65, 75, 85, 100], 
		"heterogeneity": "Slightly", 
		"smoothness": {"age_start": 0, "amount": "Slightly", "age_end": 100}
	  }	
	  "ages": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 
		  25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 
		  48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
		  71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 
		  94, 95, 96, 97, 98, 99, 100]
	}
	
ism.py
------

.. automodule:: ism
	:members:

fit.py
------

.. automodule:: fit
	:members:

graphics.py
-----------

.. automodule:: graphics
	:members:
	:exclude-members: plot_viz_of_stochs, tally_stochs