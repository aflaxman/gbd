======================
Population Data Server
======================

The Population Data Server is necessary for DisMod, and also
interesting as a stand-alone tool.  It stores, displays, and
interpolates population estimates by region (country or GBD region).

The population data server must:

1. import data from the USEABLE_IHME_GBD_POPULATION csv file
2. interpolate data from csv file to find population by age for a given region during a given time range
3. display the population pyramid graphically


Current Implementation
----------------------

* Requirement 1, Import Data.  Implemented as a Django management command:

  .. method:: gbd.population_data_server.management.commands.load_population_csv

* Requirement 2, interpolation.  Implemented in the model method:

  .. method:: gbd.population_data_server.models.Population.gaussian_process

  Work still needed to make interpolation robust.  It would be nice to
  be able to visually compare the raw data (as it appears in the csv
  file) to the interpolated values.

* Requirement 3, display.  Implemented in the controller method:

  .. function:: gbd.population_data_server.views.population_show

  Plenty of work still needed on this.


Test Driven Development
-----------------------

This simple Django App has unit and functional tests::

    gbd.population_data_server.tests

