from django.conf.urls.defaults import *

urlpatterns = patterns('gbd.covariate_data_server.views',
    (r'summary/(\d+)$', 'covariate_type_show'),
    (r'show/([\w-]+)\+(\w+)\+(\w+).(\w+)$', 'covariate_show'),
    (r'upload/$', 'covariate_upload'),
)
