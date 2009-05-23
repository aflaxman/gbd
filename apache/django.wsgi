import os, sys

sys.path.append('/usr/local/lib/python2.5/site-packages')
sys.path.append('/net/gs/vol1/home/abie/lib/simplejson-2.0.3-py2.5-linux-x86_64.egg')
sys.path.append('/net/gs/vol1/home/abie/lib/twill-0.9-py2.4.egg')
sys.path.append('/home/abie')
sys.path.append('/home/abie/gbd')

os.environ['DJANGO_SETTINGS_MODULE'] = 'gbd.settings'

import django.core.handlers.wsgi

application = django.core.handlers.wsgi.WSGIHandler()


