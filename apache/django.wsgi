import os
import sys

os.environ['DJANGO_SETTINGS_MODULE'] = 'gbd.settings'

sys.path.append('/home/abie/gbd')
import django.core.handlers.wsgi
application = django.core.handlers.wsgi.WSGIHandler()