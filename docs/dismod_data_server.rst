==================
DisMod Data Server
==================

The DisMod Data Server will be responsible for all of the web-based
generic disease modeling.  It will have four major functions: Loading
Data, Viewing Models, Running Models, and Adjusting Models.
Requirements for each of these functions are detailed below, as well
as additional requirements.


Load Data
---------

1.  Accept csv in `Standard file format <file_formats.html>`_

2.  Provide an informative error message if the csv is not parsable,
    including the row number where the error appears (for easy correction).

3.  Generate Disease Model with all data included, adding country
    level covariates and age weights, when necessary and available.

4.  Load Regional Mortality Curves from a csv file, and easily merge
    them into Disease Models.

5.  Ability to add additional data to an existing model.

View Model
----------

1.  Display data (as plots and as tables) as well as priors and estimates by:

    * type
    * year
    * sex
    * region

2.  Organize display as:

    * Panels
    * Overlay
    * Sparklines with zoom
    * Stack
    * Geodata

3.  Web-based, and on reload show the most recent version of currently
    shown condition, [region], [year], [sex].

4.  Version history with simple navigation through previous models (to
    see different estimates for different priors, etc.)

5.  A way to select and inspect individual data points.


Run Model
---------

1.  Can be run locally (from ipython shell) or on the IHME cluster,
    via web interface.

2.  Can run on a subset of the data quickly (for exploratory
    development of priors). Could be implemented as "run all
    region/year/sexes individually", "borrow strength within region",
    "borrow strength across regions".

3.  Can produce a `json description of the model <dismod_data_json.html>`_ in
    response to an HTTP GET request.

4.  Can update a model in the database based on a `json description of
    the model <dismod_data_json.html>`_ in an HTTP POST request.

5.  Has the ability to "Warm Start", i.e. use results of previous
    calculation as initial values to speed convergence.

6.  Ability to produce json list of ids for all disease models with
    needs_to_run flag set to True, in response to a get request to
    ``/dismod/job_queue/list/?format=json``.

7.  Ability to set disease model needs_to_run flags to False in
    response to a POST request to ``/dismod/job_queue/remove/<id>``.

Adjust Data and Priors
----------------------

1.  Data is adjustable, and these changes automatically track who made
    the change, when and why (by asking for the reason for the
    update).

2.  Data changes have version history so that it is easy to revert
    changes.

3.  Know what estimates must be rerun when data has changed (track
    dependency structure)

4.  Priors can be specified using:

    * text field on a webpage
    * ipython shell
    * GUI with widgets for each prior type, including:

        * smooth
        * zero
        * value
        * confidence

5.  Clear information about the meaning of each type of prior

6.  Can set priors for each data type, for each region, for each sex,
    for each year, and setting cascade in the appropriate way if they
    are not set, e.g. global priors apply to each region unless they
    are over-ridden on a region-by-region basis.

Export
------

1.  Export plots, spreadsheets, and maps (to be filled in later to match what we've already done)

2.  Export country-level samples from the posterior distribution, for a specified parameter type, age range (with specified age-standardization), year, and sex.  (This could be a python command only, i.e. not in GUI)

Implementation
--------------

Loading data is implemented, and should be working, but does not have
covariates yet, because the covariate data server still needs to be
written.

Viewing data is implemented, and works acceptible.  Speed increase is
possible, by caching images on disk, instead of re-rendering for every
page view.  Visual improvements, like line width and alpha levels may
still be needed.

Running data is implemented via the command line in ipython, and via
the web.

Setting priors is has a preliminary web interface.

Test Driven Development
-----------------------

There should be unit and functional tests to cover all the code in
this module.

Loading data has decent test coverage, although the informative error
messaging is not covered very well.  There is a failing test for a
covariate to remind me where this funcion is missing.

Viewing and running have pretty good test coverage.

Prior setting has basic test coverage, but does not test how errors in
form submission are handled.