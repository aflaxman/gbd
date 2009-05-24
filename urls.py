from django.conf.urls.defaults import *

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    (r'^admin/(.*)', admin.site.root),
    (r'^accounts/login/$', 'django.contrib.auth.views.login'),
    (r'^$', 'gbd.views.index'),
                       
    (r'^public/$', 'django.views.static.serve',
        {'document_root': 'public',
         'path': 'index.html'}),
    (r'^public/(?P<path>.*)$', 'django.views.static.serve',
        {'document_root': 'public'}),

    (r'^population/', include('gbd.population_data_server.urls')),

    (r'^new/', include('gbd.new_dm3.urls')),
    (r'^', include('gbd.dismod3.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    (r'^admin/doc/', include('django.contrib.admindocs.urls')),

)
