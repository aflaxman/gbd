from django.contrib import admin

from gbd.population_data_server.models import *

admin.site.register(Population, PopulationAdmin)

