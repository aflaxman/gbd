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
    (r'show/overlay_(\d+)_([\w-]+)\+([\w-]+)\+([\w-]+)\+(\w+)\+(\w+)\.(\w+)', 'dismod_plot', {'style': 'overlay'}),
    (r'show/bar_(\d+)_([\w-]+)\+([\w-]+)\+([\w-]+)\+(\w+)\+(\w+)\.(\w+)', 'dismod_plot', {'style': 'bar'}),
    (r'show/sparkline_(\d+)_([\w-]+)\+([\w-]+)\+([\w-]+)\+(\w+)\+(\w+)\.(\w+)', 'dismod_plot', {'style': 'sparkline'}),
    (r'show/tile_(\d+)_([\w-]+)\+([\w-]+)\+([\w-]+)\+(\w+)\+(\w+)\.(\w+)', 'dismod_plot', {'style': 'tile'}),
    (r'show/map_(\d+)$', 'dismod_show_map'),
    (r'show/(\d+)$', 'dismod_show'),

    (r'show/(\d+)/([\w-]+)$', 'dismod_show_by_region'),
    (r'show/(\d+)/([\w-]+)\.(\w+)$', 'dismod_show_by_region'),
    
    (r'show/plot_selected_regions_(\d+)$', 'dismod_show_selected_regions'),
    (r'show/plot_all_years_(\d+)$', 'dismod_show_all_years'),
    (r'show/plot_all_sexes_(\d+)$', 'dismod_show_all_sexes'),

    (r'show/(\d+)\.(\w+)$', 'dismod_show'),
    (r'show/([\w-]+)$', 'dismod_find_and_show'),
    (r'show/([\w-]+)\.(\w+)$', 'dismod_find_and_show'),
    (r'show/$', 'dismod_list'),

    (r'compare/$', 'dismod_compare'),
    (r'compare/comparison_(\d+)_(\d+)_([\w\+-]+).(\w+)$', 'dismod_comparison_plot'),

    (r'summary/(\d+)$', 'dismod_summary'),
    (r'export/(\d+)$', 'dismod_export'),

    (r'emp_priors/(\d+)$', 'dismod_show_emp_priors'),
    (r'emp_priors/(\d+).csv$', 'dismod_show_emp_priors', {'format': 'csv'}),
    (r'emp_priors/alpha_(\d+).(\w+)$', 'dismod_show_emp_priors', {'effect': 'alpha'}),
    (r'emp_priors/beta_(\d+).(\w+)$', 'dismod_show_emp_priors', {'effect': 'beta'}),
    (r'emp_priors/gamma_(\d+).(\w+)$', 'dismod_show_emp_priors', {'effect': 'gamma'}),
    (r'emp_priors/delta_(\d+).(\w+)$', 'dismod_show_emp_priors', {'effect': 'delta'}),

    (r'update_covariates/(\d+)', 'dismod_update_covariates'),
    (r'set_covariates/(\d+)', 'dismod_set_covariates'),

    (r'adjust_priors/(\d+)', 'dismod_adjust_priors'),
    (r'preview_prior/(\d+)', 'dismod_preview_priors'),
    (r'run/(\d+)', 'dismod_run'),
    (r'show_status/(\d+)$', 'dismod_show_status'),
    (r'init_log/(\d+)/([\w-]+)/(\d+)$', 'dismod_init_log'),
    (r'log_status/(\d+)/([\w-]+)/([\w-]+)/(\w+)$', 'dismod_log_status'),
    (r'server_load$', 'dismod_server_load'),
    )
