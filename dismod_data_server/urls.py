from django.conf.urls.defaults import *

urlpatterns = patterns(
    'gbd.dismod_data_server.views',

    (r'data/upload/$', 'data_upload'),
    (r'data/upload/(\d+)$', 'data_upload'),
    (r'data/(\d+)$', 'data_show'),
    (r'data/(\d+).(\w+)$', 'data_show'),
                       
    (r'upload/$', 'dismod_upload'),
    (r'job_queue/list/', 'job_queue_list'),
    (r'job_queue/add/(\d+)', 'job_queue_add'),
    (r'job_queue/remove', 'job_queue_remove'),
    
    (r'show/spark_(\d+)\.(\w+)$', 'dismod_sparkplot'),
    (r'show/overlay_(\d+)_([\w-]+)\+([\w-]+)\+([\w-]+)\+(\w+)\+(\w+)\.(\w+)', 'dismod_overlay_plot'),
    (r'show/tile_(\d+)_([\w-]+)\+([\w-]+)\+([\w-]+)\+(\w+)\+(\w+)\.(\w+)', 'dismod_tile_plot'),
    (r'show/(\d+)$', 'dismod_show'),
    (r'show/(\d+)/([\w-]+)/(\w+)/(\w+)$', 'dismod_show_by_region_year_sex'),
    (r'show/(\d+)/([\w-]+)/(\w+)$', 'dismod_show_by_region'),
    (r'show/(\d+)/([\w-]+)/(\w+)/(\w+)/(\w+)$', 'dismod_show_by_region_year_sex'),
    (r'show/(\d+)/([\w-]+)$', 'dismod_show_by_region'),
    (r'show/(\d+)\.(\w+)$', 'dismod_show'),
    (r'show/([\w-]+)$', 'dismod_find_and_show'),
    (r'show/([\w-]+)\.(\w+)$', 'dismod_find_and_show'),
    (r'show/$', 'dismod_list'),

    (r'summary/(\d+)$', 'dismod_summary'),
    (r'export/(\d+)$', 'dismod_export'),

    (r'update_covariates/(\d+)', 'dismod_update_covariates'),
    (r'adjust/(\d+)', 'dismod_adjust'),
    (r'preview_prior/(\d+)', 'dismod_preview_priors'),
    (r'run/(\d+)', 'dismod_run'),
    )
