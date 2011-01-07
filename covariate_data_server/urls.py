from django.conf.urls.defaults import *

urlpatterns = patterns('gbd.covariate_data_server.views',
    (r'summary/(\d+)$', 'covariate_type_show'),
    (r'show/([\w-]+)\+(\w+)\+(\w+).(\w+)$', 'covariate_show'),
    (r'upload/$', 'covariate_upload'),
    (r'type/$', 'covariate_type_list_show'),
    (r'datacount/(\d+)$', 'covariate_data_count_show'),
    (r'show/([\w-]+)\+(\w+).(\w+)$', 'covariate_data_value_show'),
    (r'notes/(\d+)$', 'covariate_notes_show'),
)
