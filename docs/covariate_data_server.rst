
=====================
Covariate Data Server
=====================

The Covariate Data Server is similar to the Population Data Server,
but has no age-specific information.  It stores and displays estimates
of relevant covariates by country, year, and optionally sex (e.g. GDP or education).

The covariate data server must:

1. import data from csv files, with columns 'iso3', 'year', 'sex' (optional), ‘age’ (optional) and <covariate name>

2. record metadata about the upload: who, when, where the data came from, when it was last modified, and notes including a description, how the data was created, and how the missing values were filled in

3. show the amount of data by country, highlighting the countries without all years of data

4. show the data itself, for visual inspection

5. serve the data to other dismod components, for example when
   importing dismod disease data

6. serve simple transformations of the data: log, logit, squared, cubed, lag/lead by n years, normalized, quantized


Current Implementation
----------------------

There are now some stub views and command-line tools.  To upload a new
country-level covariate, prepare a csv file for all country-years with
columns as listed above, and then use the url http://winthrop.ihme.washington.edu/covariate/upload

This view requires you input the name of the column containing the covariate of interest, for example if your table looks like this:

+------+------+------+
| iso3 | year | gdp  |
+------+------+------+
| USA  | 2005 | 1.30 |
+------+------+------+
| ...  | ...  | ...  |
+------+------+------+

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
