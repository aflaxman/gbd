from django.conf.urls.defaults import *

urlpatterns = patterns(
    'gbd.dismod_data_server.views',

    (r'data/upload/$', 'data_upload'),
    (r'data/upload/(\d+)$', 'data_upload'),
    (r'data/(\d+)$', 'data_show'),
    (r'data/(\d+).(\w+)$', 'data_show'),
                       
    (r'upload/$', 'dismod_upload'),
    (r'show/spark_(\d+)\.(\w+)$', 'dismod_sparkplot'),
    (r'show/overlay_(\d+)_([\w-]+)\+([\w-]+)\+([\w-]+)\+(\w+)\+(\w+)\.(\w+)', 'dismod_overlay_plot'),
    (r'show/tile_(\d+)_([\w-]+)\+([\w-]+)\+([\w-]+)\+(\w+)\+(\w+)\.(\w+)', 'dismod_tile_plot'),
    (r'show/(\d+)$', 'dismod_show'),
    (r'show/(\d+)\.(\w+)$', 'dismod_show'),
    (r'show/(\w+)$', 'dismod_find_and_show'),
    (r'show/$', 'dismod_list'),
    )
