from gbd.dismod3.models import *
from django.contrib import admin

admin.site.register(Disease)
admin.site.register(Region)
#admin.site.register(Study)
#admin.site.register(Table)
admin.site.register(Rate, RateAdmin)
admin.site.register(Population, PopulationAdmin)
admin.site.register(AgeSpecificRateFunction, ASRFAdmin)
admin.site.register(DiseaseModel)
#admin.site.register(Jobs)
