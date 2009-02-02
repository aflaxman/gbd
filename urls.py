from django.conf.urls.defaults import *

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    (r'^admin/(.*)', admin.site.root),
    (r'^public/(?P<path>.*)$', 'django.views.static.serve',
        {'document_root': '/home/abie/dismod/gbd/public'}),
    (r'^public/(?P<path>.*)$', 'django.views.static.serve',
        {'document_root': '/home/a/repo2/dismod/dismod3/public'}),
    (r'^', include('gbd.dismod3.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

)
