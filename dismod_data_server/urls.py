from django.conf.urls.defaults import *

urlpatterns = patterns('gbd.dismod_data_server.views',
    (r'show/(\d+)$', 'dismod_show'),
    (r'show/(\d+).(\w+)$', 'dismod_show'),

    (r'data/upload/$', 'data_upload'),
    (r'data/(\d+)$', 'data_show'),
    (r'data/(\d+).(\w+)$', 'data_show'),
)
