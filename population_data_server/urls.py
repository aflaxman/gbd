from django.conf.urls.defaults import *

urlpatterns = patterns('gbd.population_data_server.views',
    (r'show/(\d+)$', 'population_show'),
)
