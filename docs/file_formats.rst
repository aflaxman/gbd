====================================================
Proposed GBD2005 Data Submission File Specifications
====================================================

-------
Purpose
-------

The GBD2005 project aims to estimate the disease burden
associated with more than 200 diseases and their disabling sequelae as
well as the size of disease burden that can be attributed to major
risk factors.  Expert groups are currently collecting as much
information as possible on: (a) the occurrence of diseases and
sequelae (prevalence and incidence) and other epidemiological disease
parameters (remission, mortality risk, average duration, severity
distribution); (b) prevalence of exposure to risk factors; and (c) the
risk of disease by level of exposure.  From the 'raw' data collected
during the systematic reviews, estimates will be derived for 21 World
regions by age and sex with uncertainty intervals.  This requires a
number of imputations to deal with missing values and to check
available estimates for internal consistency.  An updated version of
DisMod is key in this process for disease estimation.  DisMod will
have a 'data pre-processing' function to: (a) determine common age
patterns; and (b) impute regional estimates from the data supplied by
the Expert Groups from their systematic reviews.  Next, it will allow
checks of internal consistency between the available disease
parameters and, finally, it will produce an internally-consistent full set of
epidemiological parameters describing a disease and its disabling
sequelae.

The Comparative Risk Assessment will use similar imputation techniques
to derive the most plausible regional estimates of risk factor
exposure; that will be outlined in a similar but separate data
submission file specification.

*Although not all diseases will use DisMod, we still propose this
format for all submissions of data.*

----------------------------
Minimal Required Information
----------------------------

For DisMod and other analysis purposes, we will need certain minimal
information about the parameter of interest (e.g. incidence,
remission, case fatality, etc) or proportion (e.g. prevalence).
Currently these are: GBD Cause, GBD Cause Code, Sequela, Case
Definition, Region, Parameter, Sex, Country, Urbanicity, Coverage, Age
Start, Age End, Age Units, Estimate Year Start, Estimate Year End,
Parameter Value, Standard Error, Units, Type of Confidence Limit,
Sampling Strategy, Total Study Size N, Study Information, and
Citation.

----------------------
Additional Information
----------------------

There is disease-specific additional information that expert groups
will want to include in their models, which may be used to determine
the plausible estimates for GBD Regions without sufficient study data
available or to control for excess variance. This information may vary
depending on the particular parameter or disease in question or
describe the population from which the information is derived, but
some examples are: study type (cross-sectional vs cohort vs ...),
measurement technique (e.g. various diagnostics for diabetes), and
quality of study (representativeness, case definition, confounding
etc.).

Most disease parameters will have an age pattern that is consistent
between studies from the same region or globally. Obviously, rates
from different age ranges in the same study would suggest an age
pattern much more strongly than if they come from different studies.
Therefore, ideally, we want information on each parameter in the
greatest detail possible by age and sex.

--------------------
Proposed File Format
--------------------

For each cause, the expert group collects the parameter information as
rows in a csv file (or a format that can be converted to a csv file,
such as Access or Excel) with certain required columns and any
additional columns that they think might be important.  When in doubt,
the advice is to **leave variables *in* as additional optional columns**.
It is easier to ignore some variables than to add additional variables
later.


Required data fields:

+---------------------------------+--------+------------------------------+
|Name                             |  Type  |  Limit                       |
+---------------------------------+--------+------------------------------+
|GBD Cause                        |  str   |  one of the GBD causes       |
+---------------------------------+--------+------------------------------+
|Region                           |  str   |  one of the GBD regions      |
+---------------------------------+--------+------------------------------+
|Parameter                        |  str   |  standardize_data_type       |
+---------------------------------+--------+------------------------------+
|Sex                              |  str   |  standardize_sex             |
+---------------------------------+--------+------------------------------+
|Country ISO3 Code                |  str   |  an ISO3 code in the region  |
+---------------------------------+--------+------------------------------+
|Age Start                        |  int   |  [0, 100], <= Age End        |
+---------------------------------+--------+------------------------------+
|Age End                          |  int   |  [0, 100], >= Age Start      |
+---------------------------------+--------+------------------------------+
|Year Start                       |  int   |  [1980, 2010], <= Year End   |
+---------------------------------+--------+------------------------------+
|Year End                         |  int   |  [1980, 2010], >= Year Start |
+---------------------------------+--------+------------------------------+
|Parameter Value                  |  float |  >= 0                        |
+---------------------------------+--------+------------------------------+
|Units                            |  float |  >= 1                        |
+---------------------------------+--------+------------------------------+

Recommended data fields:

+---------------------------------+--------+------------------------------+
|Name                             |  Type  |  Limit                       |
+---------------------------------+--------+------------------------------+
|Study ID                         |  int   |  >= 0                        |
+---------------------------------+--------+------------------------------+
|Sequela                          |  str   |  one of the GBD sequela codes|
+---------------------------------+--------+------------------------------+
|Case Definition                  |  str   |  none                        |
+---------------------------------+--------+------------------------------+
|Coverage                         |  float |  [0,1]                       |
+---------------------------------+--------+------------------------------+
|Study Size N For This Year & Sex |  int   |  > 0, <= Total Study Size N  |
+---------------------------------+--------+------------------------------+
|Lower CI                         |  float |  > 0 <= Parameter Value      |
+---------------------------------+--------+------------------------------+
|Upper CI                         |  float |  >= Parameter Value          |
+---------------------------------+--------+------------------------------+
|Standard Error                   |  float |  > 0                         |
+---------------------------------+--------+------------------------------+
|Total Study Size N               |  int   |  > 0                         |
+---------------------------------+--------+------------------------------+
|Design Factor                    |  float |  >= 1                        |
+---------------------------------+--------+------------------------------+
|Citation                         |  str   |  none                        |
+---------------------------------+--------+------------------------------+
|Urbanicity                       |  float |  [0, 1]                      |
+---------------------------------+--------+------------------------------+


Optional data fields: No checks         
