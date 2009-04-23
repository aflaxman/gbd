from django.conf.urls.defaults import *

urlpatterns = patterns('new_dm3.views',
    (r'data/(\d+)$', 'data_show'),
    (r'data/new$', 'data_new'),

    (r'dm/(\d+)\.(\w+)$', 'disease_model_show'),
    (r'dm/(\d+)$', 'disease_model_show'),
)
