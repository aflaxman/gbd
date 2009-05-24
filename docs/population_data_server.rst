
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
======================

Requirement 1
-------------
Implemented as a Django management command:

.. automodule:: gbd.population_data_server.management.commands.load_population_csv

Requirement 2
-------------
Implemented in the model method:

.. automethod:: gbd.population_data_server.models.Population.gaussian_process

Requirement 3
-------------
Implemented in the controller method:

.. function:: gbd.population_data_server.views.population_show


Test Driven Development
=======================

This simple Django App has unit and functional tests:

.. automodule:: gbd.population_data_server.tests
   :members:
   :undoc-members:

