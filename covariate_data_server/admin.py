from django.contrib import admin

from gbd.covariate_data_server.models import *

admin.site.register(Covariate, CovariateAdmin)

