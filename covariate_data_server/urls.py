from django.conf.urls.defaults import *

urlpatterns = patterns('gbd.covariate_data_server.views',
    (r'show/(\d+)$', 'covariate_show'),
    (r'show/(\d+).(\w+)$', 'covariate_show'),
)
