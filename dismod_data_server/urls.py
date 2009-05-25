from django.conf.urls.defaults import *

urlpatterns = patterns('gbd.dismod_data_server.views',

                       (r'data/upload/$', 'data_upload'),
                       (r'data/(\d+)$', 'data_show'),
                       (r'data/(\d+).(\w+)$', 'data_show'),
                       
                       (r'upload/$', 'dismod_upload'),
                       (r'show/spark_(\d+).(\w+)$', 'dismod_sparkplot'),
                       (r'show/(\d+)$', 'dismod_show'),
                       (r'show/(\d+).(\w+)$', 'dismod_show'),
                       )
