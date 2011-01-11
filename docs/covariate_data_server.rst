=====================
Covariate Data Server
=====================

The Covariate Data Server is similar to the Population Data Server,
but has no age-specific information.  It stores and displays estimates
of relevant covariates by country, year, and optionally sex (e.g. GDP or education).

The covariate data server must:

1. import data from csv files through web interface or command line interface, with columns 'iso3', 'year', 'sex' (optional), ‘age’ (optional) and <covariate name>.  instead of iso3, allow user to specify gbd_region as a column, and return this value for all countries in the region.

2. record metadata about the upload: who, when, where the data came from, when it was last modified, and notes including a description, how the data was created, and how the missing values were filled in

3. show the amount of data by country (for all sexes together), highlighting the countries without all years of data (what does it mean to have "all" years of data? should be set for the project year_start (e.g. 1950) and year_end (e.g. 2010) for covariates.)

4. show the data itself, for visual inspection

5. serve the data to other dismod components, for example when
   importing dismod disease data.  e.g. in dismod_data_server/models.py the calculate_covariate method (this calculate_covariates method could be improved for speed...  it is always used to process a batch of data, and the covariate_data_server could take this into account for more efficient database access, so we need a function that is not a member of a specific Data object, which takes a list of Data objects and stores their covariate values.  Or maybe this redundant storage is not necessary, and as long as the disease_json has the covariate data when it is created, maybe duplicating the data in the mysql database is unnecessary).  

6. serve simple transformations of the data: log, logit, squared, cubed, lag/lead by n years, normalized, quantized

Integrating the covariate data server and the data checker java app
-------------------------------------------------------------------

1. The Data Checker must request a list of available covariate types from the data server

2. The Data Checker must display a list of available covariates and a set of checkboxes for allowable transformations so the user can choose which covariates with which transformations to add

3. The Data Checker must ask the covariate data server to merge the chosen (covariate, transformation)-pairs into the data, and display the results in the data checker table

4. The Data checker must mark as error (or maybe warning?) any cells which the covariate data server could not fill in


Current Implementation
----------------------

There are now some stub views and command-line tools.  To upload a new
country-level covariate, prepare a csv file for all country-years with
columns as listed above, and then use the url http://winthrop.ihme.washington.edu/covariate/upload

This view requires you input the name of the column containing the covariate of interest, for example if your table looks like this:

+------+------+------+------+
| iso3 | year | sex  | gdp  |
+------+------+------+------+
| USA  | 2005 | male | 1.30 |
+------+------+------+------+
| ...  | ...  | ...  |      |
+------+------+------+------+

The type should be `gdp`.

Note: After loading new covariates there is some work necessary to make old models recognize it.  Example::

    dm = DiseaseModel.objects.get(id=944)
    covariates, is_new = dm.params.get_or_create(key='covariates')
    covariates.json = json.dumps({'Study_level': dm.study_level_covariates(),
                                  'Country_level': dm.country_level_covariates()})
    covariates.save()


Test Driven Development
-----------------------

There should be unit and functional tests to cover all the code in
this module.
