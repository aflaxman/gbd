import os, sys

sys.path.append('/usr/local/lib/python2.5/site-packages')
sys.path.append('/usr/local/lib/python2.5/site-packages/simplejson-2.0.9-py2.5-linux-x86_64.egg')
sys.path.append('/usr/local/lib/python2.5/site-packages/twill-0.9-py2.5.egg')

sys.path.append('/home/abie/dismod_3_beta')
sys.path.append('/home/abie/dismod_3_beta/gbd')


# set home to a directory where user apache can created .matplotlib dir
os.environ['HOME'] = '/usr/local/apache2/logs'

# tell django where to find the settings for the project
os.environ['DJANGO_SETTINGS_MODULE'] = 'gbd.settings'
import django.core.handlers.wsgi

application = django.core.handlers.wsgi.WSGIHandler()
