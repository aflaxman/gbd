from django.contrib import admin

from gbd.dismod_data_server.models import *

admin.site.register(Data, DataAdmin)
admin.site.register(DiseaseModel, DiseaseModelAdmin)

