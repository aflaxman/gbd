======================
Population Data Server
======================

The Population Data Server is necessary for DisMod, and also
interesting as a stand-alone tool.  It stores, displays, and
interpolates population estimates by region (country or GBD region).

The population data server must:

1. import data from the USEABLE_IHME_GBD_POPULATION csv file

2. aggregate data over countries to find the population of each GBD region

3. interpolate data from csv file to find population by age for a
   given region during a given time range

4. display the population pyramid graphically


Current Implementation
----------------------

* Requirement 1, importation.  Implemented as a Django management command::

    $ python2.5 manage.py load_population_csv USABLE_IHME_GBD_POPULATION_1950-2050.csv

* Requirement 2, aggregation.  Included in the ``load_population_csv``
  management command mentioned above.

* Requirement 3, interpolation.  Implemented using PyMC Gaussian
  Processes, as a method in the ``models.Population`` model::

    >>> pop = Population.objects.latest('id')
    >>> M,C = pop.gaussian_process()
    >>> M(range(100)) # interpolated over ages [0, 1, 2, ..., 99]

  Work is still needed to make interpolation robust.  It would be nice to
  be able to visually compare the raw data (as it appears in the USABLE_IHME csv
  file) to the interpolated values.

* Requirement 4, display.  Implemented as a Django method in the
  ``views.population_show``:

  .. function:: gbd.population_data_server.views.population_show

  Plenty of work still needed on this.


Test Driven Development
-----------------------

This simple Django App has unit and functional tests::

    gbd.population_data_server.tests

