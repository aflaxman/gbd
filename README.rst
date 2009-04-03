Introduction
============

This is the currently in-development version of DisMod III, the
Generic Disease Modeling System.  It is being developed for the
Institute of Health Metrics and Evaluation at the University of
Washington, and will one day be part of the Global Burden of Disease
Study.  To learn more, visit http://www.globalburden.org/

Installation
============

Requirements: Python2.5, easy_install, numpy, scipy, matplotlib, PyMC,
Django, simplejson, sphinx, pygments, [more?]

Installation for Ubuntu::

    git clone git://github.com/aflaxman/gbd.git
    sudo apt-get install python2.5 python-setuptools
    sudo easy_install pymc django simplejson
    python2.5 manage.py syncdb

Installation may also be possible for Windows::
    http://www.enthought.com/products/epd.php
    http://code.google.com/p/pymc/downloads/list
    http://wiki.thinkhole.org/howto:django_on_windows

Configuration:
    Edit gbd/settings.py, and make all the entries of the TEMPLATE_DIRS more accurate

Testing::
    python2.5 manage.py test

Running the server (for development and interactive testing)::
    python2.5 manage.py runserver winthrop.gs.washington.edu:5432

Testing
=======

DisMod III has good test coverage.  To run the automatic tests::
    python2.5 manage.py test

There are also some tools for interactive testing, which should be
useful during model development and validation.  To enter the DisMod
shell::
    python2.5 manage.py shell
To get models in the namespace::
    from dismod3.models import *
To create a test rate function::
    import dismod3.tests.bayesian_probability_test as bpt
    rf = bpt.create_test_asrf('(age/100.)**2 / 5.')
    
Note that the shell does not reload changes automatically!  If you
change code (in bayesian_probability_test.py, say) and want to see the
effects in the shell::
    reload(bpt)
    
Migrations
==========

Changes to the model schema are somewhat difficult in Django.  Here
are some notes on how to make it a little bit easier::
    python2.5 manage.py dumpdata dismod3 >dm_data_YYYY_MM_DD.json

Make changes to the schema, for example::
    --- dismod3/models/rate.py  (revision 392)
    +++ dismod3/models/rate.py  (working copy)
    @@ -36,17 +36,29 @@
         epoch_end = models.IntegerField()
         numerator = models.IntegerField()
         denominator = models.IntegerField()
    +    rate = models.FloatField(editable=False, default=0.)
         params_json = models.TextField(default=json.dumps({}))

Drop the application tables, and then syncdb to load the migrated
tables::
    python2.5 manage.py sqlclear dismod3 |python2.5 manage.py dbshell
    python2.5 manage.py syncdb

Repopulated the database with the data you dumped::
    python2.5 manage.py loaddata dismod3 dm_data_YYYY_MM_DD.json
