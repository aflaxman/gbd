from django.core.management.base import BaseCommand, CommandError

from dismod3.models import AgeSpecificRateFunction, DiseaseModel
from dismod3.tests.bayesian_probability_test import create_test_asrf
from dismod3.bayesian_models import fit_rate_function
from dismod3.bayesian_models import fit_disease_model

class Command(BaseCommand):
    help = 'Generate rate function data and fit it with a variety of priors'

    def handle(self, *args, **options):
        if args[0] == 'dm':
            r = create_test_asrf('.05',
                                 rate_type='remission data',
                                 priors='smooth 10.0\nconfidence 1000 .1')
            p = create_test_asrf('.2+.3*.01*age',
                                 rate_type='prevalence data',
                                 priors='smooth 10.0\nincreasing 0 100\nconfidence 1000 .1')
            cf = create_test_asrf('0.',
                                  rate_type='case fatality data',
                                  priors='smooth 10.0\nconfidence 200 .1')

            dm = DiseaseModel(disease=cf.disease, region=cf.region, sex=cf.sex)
            dm.save()
            dm.rates = [r,p,cf]
            dm.save()

            print "Fitting %s\n  http://winthrop.gs.washington.edu:5432%s" % (dm.id, dm.get_absolute_url())
            fit_disease_model.mcmc_fit(dm)



            i = create_test_asrf('0.',
                                 rate_type='incidence data',
                                 priors='smooth 10.0\nconfidence 1000 .1')
            p = create_test_asrf('.5 - .3*.01*age',
                                 rate_type='prevalence data',
                                 priors='smooth 10.0\ndecreasing 0 100\nconfidence 1000 .1')
            cf = create_test_asrf('0.',
                                  rate_type='case fatality data',
                                  priors='smooth 10.0\nconfidence 1000 .1')

            dm = DiseaseModel(disease=cf.disease, region=cf.region, sex=cf.sex)
            dm.save()
            dm.rates = [i,p,cf]
            dm.save()

            print "Fitting %s\n  http://winthrop.gs.washington.edu:5432%s" % (dm.id, dm.get_absolute_url())
            fit_disease_model.mcmc_fit(dm)


            i = create_test_asrf('0.',
                                 rate_type='incidence data',
                                 priors='smooth 10.0\nconfidence 1000 .1')
            r = create_test_asrf('.05',
                                 rate_type='remission data',
                                 priors='smooth 10.0\nconfidence 1000 .1')
            cf = create_test_asrf('0.',
                                  rate_type='case fatality data',
                                  priors='smooth 10.0\nconfidence 1000 .1')

            dm = DiseaseModel(disease=cf.disease, region=cf.region, sex=cf.sex)
            dm.save()
            dm.rates = [i,r,cf]
            dm.save()

            print "Fitting %s\n  http://winthrop.gs.washington.edu:5432%s" % (dm.id, dm.get_absolute_url())
            fit_disease_model.mcmc_fit(dm)
            i = create_test_asrf('.1',
                                 rate_type='incidence data',
                                 priors='smooth 10.0\nconfidence 1000 .1')
            r = create_test_asrf('.05',
                                 rate_type='remission data',
                                 priors='smooth 10.0\nconfidence 1000 .1')
            p = create_test_asrf('.07*(age/100.0)**2',
                                 rate_type='prevalence data',
                                 priors='smooth 10.0\nconfidence 1000 .1')

            dm = DiseaseModel(disease=cf.disease, region=cf.region, sex=cf.sex)
            dm.save()
            dm.rates = [i,r,p]
            dm.save()

            print "Fitting %s\n  http://winthrop.gs.washington.edu:5432%s" % (dm.id, dm.get_absolute_url())
            fit_disease_model.mcmc_fit(dm)

            print "\nModel %d is fit\n  http://winthrop.gs.washington.edu:5432%s" % (dm.id, dm.get_absolute_url())
        elif args[0] == 'rf':
            if args[1] == 'unimodal':
                rf = create_test_asrf('.6-2.*(.5-age/100.0)**2',
                                     rate_type='incidence data',
                                     priors='smooth 10.0\nunimodal 0 100\nconfidence 1000 .1')
                print "Fitting %s\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())
                fit_rate_function.mcmc_fit(rf)
                print "\nModel %d is fit\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())

                rf = create_test_asrf('.06-.2*(.5-age/100.0)**2',
                                     rate_type='incidence data',
                                     priors='smooth 10.0\nunimodal 0 100')
                print "Fitting %s\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())
                fit_rate_function.mcmc_fit(rf)
                print "\nModel %d is fit\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())
            elif args[1] == 'increasing':
                rf = create_test_asrf('.01 + .5*(age/100.0)**2',
                                     rate_type='remission data',
                                     priors='smooth 10.0\nincreasing 0 100\nconfidence 1000 .1')
                print "Fitting %s\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())
                fit_rate_function.mcmc_fit(rf)
                print "\nModel %d is fit\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())
                
                rf = create_test_asrf('.001 + .05*(age/100.0)',
                                     rate_type='incidence data',
                                     priors='smooth 10.0\nincreasing 0 100')
                print "Fitting %s\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())
                fit_rate_function.mcmc_fit(rf)
                print "\nModel %d is fit\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())
            elif args[1] == 'decreasing':
                rf = create_test_asrf('.61 - .5*(age/100.0)**2',
                                     rate_type='remission data',
                                     priors='smooth 10.0\ndecreasing 0 100')
                print "Fitting %s\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())
                fit_rate_function.mcmc_fit(rf)
                print "\nModel %d is fit\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())
                
                rf = create_test_asrf('.061 - .05*(age/100.0)',
                                     rate_type='incidence data',
                                     priors='smooth 10.0\ndecreasing 0 100')
                print "Fitting %s\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())
                fit_rate_function.mcmc_fit(rf)
                print "\nModel %d is fit\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())

