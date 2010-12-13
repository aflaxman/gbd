#!/usr/bin/python2.5
""" Distribute the parameter estimation throughout our cluster.

gbd_fit has options to start in daemon mode or in fitting mode.

Examples
--------

$ python gbd_fit --daemon    # launch daemon that will fit models as they become available
$ python gbd_fit 10   # launch fitting calculation to estimate parameters for model #10
$ python gbd_fit 10 --nofit -t incidence -p 'smooth 25'  # set the hyper-prior on incidence to 'smooth 25' and save it, without running the model

"""

import time
import optparse
import subprocess
import signal

import dismod3
from dismod3.utils import clean, gbd_keys, type_region_year_sex_from_key
from dismod3.plotting import GBDDataHash

import sys
from os import popen
import os
from shutil import rmtree
import daemon

def tweet(message,
          user=dismod3.settings.DISMOD_TWITTER_NAME,
          password=dismod3.settings.DISMOD_TWITTER_PASSWORD):
    print 'tweeting %s for %s' % (message, user)

    message = '#dismod %s' % message
 
    url = 'http://twitter.com/statuses/update.xml' 
    curl = '#curl -s -u %s:%s -d status="%s" %s' % (user,password,message,url)
    try:
        pipe = popen(curl, 'r')
    except:
        pass

def log(message):
    sys.stdout.write('%s  gbd_fit  %s\n' % (time.strftime("%Y-%m-%d %H:%M:%S"), message))
    sys.stdout.flush()

def term(self, *args):
    log('dismod3 daemon received SIGTERM')
    sys.exit()

def main():
    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--type', dest='type',
                      help='only estimate given parameter type (valid settings ``incidence``, ``prevalence``, ``remission``, ``excess-mortality``) (emp prior fit only)')
    parser.add_option('-s', '--sex', dest='sex',
                      help='only estimate given sex (valid settings ``male``, ``female``, ``all``)')
    parser.add_option('-y', '--year',
                      help='only estimate given year (valid settings ``1990``, ``2005``)')
    parser.add_option('-r', '--region',
                      help='only estimate given GBD Region')

    parser.add_option('-d', '--daemon',
                      action='store_true', dest='daemon')

    parser.add_option('-l', '--log',
                      action='store_true', dest='log',
                      help='log the job running status')

    (options, args) = parser.parse_args()

    if options.daemon:
        if len(args) != 0:
            parser.error('incorrect number of arguments for daemon mode (should be none)')

        #daemon.daemonize('/dev/null', dismod3.settings.DAEMON_LOG_FILE, dismod3.settings.DAEMON_LOG_FILE)
        f = open(dismod3.settings.GBD_FIT_LOCK_FILE, 'w')
        f.write(str(os.getpid()))
        f.close()
        signal.signal(signal.SIGTERM, term)
        try:
            #tweet('starting dismod3 daemon...')
            log('starting dismod3 daemon...')
            daemon_loop()
        finally:
            #tweet('...dismod3 daemon shutting down')
            log('dismod3 daemon shutting down')

    else:
        if len(args) != 1:
            parser.error('incorrect number of arguments')

        try:
            id = int(args[0])
        except ValueError:
            parser.error('disease_model_id must be an integer')

        fit(id, options)

def daemon_loop():
    on_sge = dismod3.settings.ON_SGE
    while True:
        try:
            job_queue = dismod3.get_job_queue()
        except:
            job_queue = []
        
        for param_id in job_queue:
            #tweet('processing job %d' % id)
            log('processing job %d' % param_id)
            job_params = dismod3.remove_from_job_queue(param_id)
            id = int(job_params['dm_id'])
            dm = dismod3.get_disease_model(id)

            # make a working directory for the id
            dir = dismod3.settings.JOB_WORKING_DIR % id
            if os.path.exists(dir):
                dismod3.disease_json.random_rename(dir)
            os.makedirs(dir)

            estimate_type = dm.params.get('run_status', {}).get('estimate_type', 'fit all individually')

            # sort the regions so that the data rich regions are fit first
            #data_hash = GBDDataHash(dm.data)
            #sorted_regions = sorted(dismod3.gbd_regions, reverse=True,
                                    #key=lambda r: len(data_hash.get(region=r)))

            if estimate_type == 'Fit continuous single parameter model':
                #dismod3.disease_json.create_disease_model_dir(id)
                o = '%s/continuous_spm.stdout' % dir
                e = '%s/continuous_spm.stderr' % dir
                if on_sge:
                    print o
                    print e
                    call_str = 'qsub -cwd -o %s -e %s ' % (o, e) \
                               + 'run_on_cluster.sh /home/OUTPOST/abie/gbd_dev/gbd/fit_continuous_spm.py %d' % id
                else:
                    call_str = 'python -u /home/abie/gbd/fit_continuous_spm.py %d 2>%s |tee %s' % (id, e, o)
                subprocess.call(call_str, shell=True)
                continue
            
            if estimate_type.find('posterior') != -1:
                #fit each region/year/sex individually for this model
                regions_to_fit = dm.params.get('run_status', {}).get('regions_to_fit', [])
                if regions_to_fit[0] == 'all_regions':
                    regions_to_fit = dismod3.gbd_regions
                d = '%s/posterior' % dir
                if os.path.exists(d):
                    rmtree(d)
                os.mkdir(d)
                os.mkdir('%s/stdout' % d)
                os.mkdir('%s/stderr' % d)
                os.mkdir('%s/pickle' % d)
                dismod3.init_job_log(id, 'posterior', param_id)
                for r in regions_to_fit:
                    for s in dismod3.gbd_sexes:
                        for y in dismod3.gbd_years:
                            # fit only one region, for the time being...
                            # TODO: make region selection a user-settable option from the gui
                            #if clean(r) != 'asia_southeast':
                            #    continue
                            k = '%s+%s+%s' % (clean(r), s, y)
                            o = '%s/stdout/%s' % (d, k)
                            e = '%s/stderr/%s' % (d, k)
                            if on_sge:
                                call_str = dismod3.settings.GBD_FIT_STR % (o, e, '-l -r %s -s %s -y %s' % (clean(r), s, y), id)
                                subprocess.call(call_str, shell=True)
                            else:
                                call_str = dismod3.settings.GBD_FIT_STR % ('-l -r %s -s %s -y %s' % (clean(r), s, y), id, o, e)
                                subprocess.call(call_str, shell=True)
                            #time.sleep(1.)

            elif estimate_type.find('empirical priors') != -1:
                # fit empirical priors (by pooling data from all regions
                d = '%s/empirical_priors' % dir
                if os.path.exists(d):
                    rmtree(d)
                os.mkdir(d)
                os.mkdir('%s/stdout' % d)
                os.mkdir('%s/stderr' % d)
                os.mkdir('%s/pickle' % d)
                dismod3.init_job_log(id, 'empirical_priors', param_id)
                for t in ['excess-mortality', 'remission', 'incidence', 'prevalence']:
                    o = '%s/stdout/%s' % (d, t)
                    e = '%s/stderr/%s' % (d, t)
                    if on_sge:
                        subprocess.call(dismod3.settings.GBD_FIT_STR % (o, e, '-l -t %s' % t, id), shell=True)
                    else:
                        subprocess.call(dismod3.settings.GBD_FIT_STR % ('-l -t %s' % t, id, o, e), shell=True)

            else:
                #tweet('unrecognized estimate type: %s' % estimate_type)
                log('unrecognized estimate type: %s' % estimate_type)
            
        time.sleep(dismod3.settings.SLEEP_SECS)
        
def fit(id, opts):
    fit_str = '(%d) %s %s %s' % (id, opts.region or '', opts.sex or '', opts.year or '')
    #tweet('fitting disease model %s' % fit_str)
    sys.stdout.flush()
    
    # update job status file
    if opts.log:
        if opts.type and not (opts.region and opts.sex and opts.year):
            dismod3.log_job_status(id, 'empirical_priors', opts.type, 'Running')
        elif opts.region and opts.sex and opts.year and not opts.type:
            dismod3.log_job_status(id, 'posterior', '%s--%s--%s' % (opts.region, opts.sex, opts.year), 'Running')

    dm = dismod3.get_disease_model(id)
    fit_str = '%s %s' % (dm.params['condition'], fit_str)

    sex_list = opts.sex and [ opts.sex ] or dismod3.gbd_sexes
    year_list = opts.year and [ opts.year ] or dismod3.gbd_years
    region_list = opts.region and [ opts.region ] or dismod3.gbd_regions
    keys = gbd_keys(region_list=region_list, year_list=year_list, sex_list=sex_list)

    # fit empirical priors, if type is specified
    if opts.type:
        fit_str += ' emp prior for %s' % opts.type
        #print 'beginning ', fit_str
        import dismod3.neg_binom_model as model

        dir = dismod3.settings.JOB_WORKING_DIR % id
        model.fit_emp_prior(dm, opts.type, dbname='%s/empirical_priors/pickle/dm-%d-emp_prior-%s.pickle' % (dir, id, opts.type))

    # if type is not specified, find consistient fit of all parameters
    else:
        import dismod3.gbd_disease_model as model

        # get the all-cause mortality data, and merge it into the model
        mort = dismod3.get_disease_model('all-cause_mortality')
        dm.data += mort.data

        # fit individually, if sex, year, and region are specified
        if opts.sex and opts.year and opts.region:
            dm.params['estimate_type'] = 'fit individually'

        # fit the model
        #print 'beginning ', fit_str
        dir = dismod3.settings.JOB_WORKING_DIR % id
        model.fit(dm, method='map', keys=keys, verbose=1)
        model.fit(dm, method='mcmc', keys=keys, iter=10000, thin=5, burn=5000, verbose=1,
                  dbname='%s/posterior/pickle/dm-%d-posterior-%s-%s-%s.pickle' % (dir, id, opts.region, opts.sex, opts.year))
        #model.fit(dm, method='mcmc', keys=keys, iter=1, thin=1, burn=0, verbose=1)

    # remove all keys that have not been changed by running this model
    for k in dm.params.keys():
        if type(dm.params[k]) == dict:
            for j in dm.params[k].keys():
                if not j in keys:
                    dm.params[k].pop(j)

    # post results to dismod_data_server
    # "dumb" error handling, in case post fails (try: except: sleep random time, try again, stop after 4 tries)
    from twill.errors import TwillAssertionError
    from urllib2 import URLError
    import random

    PossibleExceptions = [TwillAssertionError, URLError]
    try:
        url = dismod3.post_disease_model(dm)
    except PossibleExceptions:
        time.sleep(random.random()*30)
        try:
            url = dismod3.post_disease_model(dm)
        except PossibleExceptions:
            time.sleep(random.random()*30)
            try:
                url = dismod3.post_disease_model(dm)
            except PossibleExceptions:
                time.sleep(random.random()*30)
                url = dismod3.post_disease_model(dm)

    # form url to view results
    #if opts.sex and opts.year and opts.region:
    #    url += '/%s/%s/%s' % (opts.region, opts.year, opts.sex)
    #elif opts.region:
    #    url += '/%s' % opts.region

    # announce completion, and url to view results
    #tweet('%s fit complete %s' % (fit_str, url))
    sys.stdout.flush()

    # update job status file
    if opts.log:
        if opts.type and not (opts.region and opts.sex and opts.year):
            dismod3.log_job_status(id, 'empirical_priors', opts.type, 'Completed')
        elif opts.region and opts.sex and opts.year and not opts.type:
            dismod3.log_job_status(id, 'posterior', '%s--%s--%s' % (opts.region, opts.sex, opts.year), 'Completed')

if __name__ == '__main__':
    main()
