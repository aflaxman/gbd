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
                      help='only estimate given parameter type (valid settings ``incidence``, ``prevalence``, ``remission``, ``case-fatality``) (emp prior fit only)')
    parser.add_option('-s', '--sex', dest='sex',
                      help='only estimate given sex (valid settings ``male``, ``female``, ``all``)')
    parser.add_option('-y', '--year',
                      help='only estimate given year (valid settings ``1990``, ``2005``)')
    parser.add_option('-r', '--region',
                      help='only estimate given GBD Region')

    parser.add_option('-p', '--prior',
                      help='prior for specified parameter type (emp prior fit only)')
    parser.add_option('-P', '--prevprior',
                      help='prevalence priors')
    parser.add_option('-I', '--inciprior',
                      help='incidence priors')
    parser.add_option('-C', '--caseprior',
                      help='case-fatality priors')
    parser.add_option('-R', '--remiprior',
                      help='remission priors')

    parser.add_option('-N', '--nofit',
                      action='store_true', dest='no_fit',
                      help='do not fit the model (save priors only)')

    parser.add_option('-d', '--daemon',
                      action='store_true', dest='daemon')

    parser.add_option('-l', '--log',
                      action='store_true', dest='log',
                      help='log the job running status')

    (options, args) = parser.parse_args()

    if len(args) == 0:
        if options.daemon:
            daemon.daemonize('/dev/null', dismod3.settings.DAEMON_LOG_FILE, dismod3.settings.DAEMON_LOG_FILE)
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
            parser.error('incorrect number of arguments')
    elif len(args) == 1:
        if options.daemon:
            parser.error('incorrect number of arguments for daemon mode (should be none)')
        try:
            id = int(args[0])
        except ValueError:
            parser.error('disease_model_id must be an integer')
        fit(id, options)
    else:
        parser.error('incorrect number of arguments')

def daemon_loop():
    on_sge = dismod3.settings.ON_SGE
    while True:
        try:
            job_queue = dismod3.get_job_queue()
        except:
            job_queue = []
        
        for id in job_queue:
            #tweet('processing job %d' % id)
            log('processing job %d' % id)
            job_params = dismod3.remove_from_job_queue(id)
            id = int(job_params['dm_id'])
            dm = dismod3.get_disease_model(id)

            # make a working directory for the id
            dir = dismod3.settings.JOB_WORKING_DIR % id
            if not os.path.exists(dir):
                os.makedirs(dir)

            estimate_type = dm.params.get('run_status', {}).get('estimate_type', 'fit all individually')

            # sort the regions so that the data rich regions are fit first
            data_hash = GBDDataHash(dm.data)
            sorted_regions = sorted(dismod3.gbd_regions, reverse=True,
                                    key=lambda r: len(data_hash.get(region=r)))
            
            if estimate_type.find('posterior') != -1:
                #fit each region/year/sex individually for this model (84 processes!)
                d = '%s/posterior' % dir
                if os.path.exists(d):
                    rmtree(d)
                os.mkdir(d)
                os.mkdir('%s/stdout' % d)
                os.mkdir('%s/stderr' % d)
                dismod3.init_job_log(id, 'posterior')
                for r in sorted_regions:
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
                            time.sleep(1.)

            elif estimate_type.find('within each region') != -1:
                # fit each region individually, but borrow strength within gbd regions
                d = '%s/within_each_region' % dir
                if os.path.exists(d):
                    rmtree(d)
                os.mkdir(d)
                os.mkdir('%s/stdout' % d)
                os.mkdir('%s/stderr' % d)
                dismod3.init_job_log(id, 'within_each_region')
                for r in sorted_regions:
                    cr = clean(r)
                    o = '%s/stdout/%s' % (d, cr)
                    e = '%s/stderr/%s' % (d, cr)
                    if on_sge:
                        subprocess.call(dismod3.settings.GBD_FIT_STR % (o, e, '-l -r %s' % cr, id), shell=True)
                    else:
                        subprocess.call(dismod3.settings.GBD_FIT_STR % ('-l -r %s' % cr, id, o, e), shell=True)

            elif estimate_type.find('across all regions') != -1:
                # fit all regions, years, and sexes together
                d = '%s/across_all_regions' % dir
                if os.path.exists(d):
                    rmtree(d)
                os.mkdir(d)
                os.mkdir('%s/stdout' % d)
                os.mkdir('%s/stderr' % d)
                o = '%s/stdout/all_regions' % d
                e = '%s/stderr/all_regions' % d
                dismod3.init_job_log(id, 'across_all_regions')
                if on_sge:
                    subprocess.call(dismod3.settings.GBD_FIT_STR % (o, e, '-l', id), shell=True)
                else:
                    subprocess.call(dismod3.settings.GBD_FIT_STR % ('-l', id, o, e), shell=True)

            elif estimate_type.find('empirical priors') != -1:
                # fit empirical priors (by pooling data from all regions
                d = '%s/empirical_priors' % dir
                if os.path.exists(d):
                    rmtree(d)
                os.mkdir(d)
                os.mkdir('%s/stdout' % d)
                os.mkdir('%s/stderr' % d)
                dismod3.init_job_log(id, 'empirical_priors')
                for t in ['case-fatality', 'remission', 'incidence', 'prevalence']:
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
    if opts.log and not opts.no_fit:
        if opts.type and not (opts.region and opts.sex and opts.year):
            dismod3.log_job_status(id, 'empirical_priors', opts.type, 'Running')
        elif opts.region and opts.sex and opts.year and not opts.type:
            dismod3.log_job_status(id, 'posterior', '%s--%s--%s' % (opts.region, opts.sex, opts.year), 'Running')
        elif opts.region and not (opts.sex and opts.year and opts.type):
            dismod3.log_job_status(id, 'within_each_region', opts.region, 'Running')
        elif not (opts.region and opts.sex and opts.year and opts.type):
            dismod3.log_job_status(id, 'across_all_regions', 'all_regions', 'Running')

    dm = dismod3.get_disease_model(id)
    fit_str = '%s %s' % (dm.params['condition'], fit_str)

    sex_list = opts.sex and [ opts.sex ] or dismod3.gbd_sexes
    year_list = opts.year and [ opts.year ] or dismod3.gbd_years
    region_list = opts.region and [ opts.region ] or dismod3.gbd_regions
    keys = gbd_keys(region_list=region_list, year_list=year_list, sex_list=sex_list)
    if not opts.region:
        keys += gbd_keys(region_list=['world'], year_list=['1997'], sex_list=['total'])

    # quick way to add/replace priors from the command line
    for rate_type, priors in [ ['prevalence', opts.prevprior], ['incidence', opts.inciprior],
                               ['remission', opts.remiprior], ['case-fatality', opts.caseprior] ]:
        if priors:
            # set priors for appropriate region-year-sex submodels
            for k in keys:
                key_type = type_region_year_sex_from_key(k)[0]
                if key_type == rate_type:
                    dm.set_priors(k, priors)

    # if opts.type is specified, also set the (hyper)-priors on the empirical prior
    if opts.type and opts.prior:
        dm.set_priors(opts.type, opts.prior)
        for k in keys:
            key_type = type_region_year_sex_from_key(k)[0]
            if key_type == opts.type:
                dm.set_priors(k, opts.prior)

    # fit empirical priors, if type is specified
    if (not opts.no_fit) and opts.type:
        fit_str += ' emp prior for %s' % opts.type
        #print 'beginning ', fit_str
        #import dismod3.logit_normal_model as model
        import dismod3.neg_binom_model as model
        model.fit_emp_prior(dm, opts.type)
        
    # if type is not specified, find consistient fit of all parameters
    elif not opts.no_fit:
        import dismod3.gbd_disease_model as model

        # get the all-cause mortality data, and merge it into the model
        mort = dismod3.get_disease_model('all-cause_mortality')
        dm.data += mort.data

        # fit individually, if sex, year, and region are specified
        if opts.sex and opts.year and opts.region:
            dm.params['estimate_type'] = 'fit individually'

        # fit the model
        #print 'beginning ', fit_str
        model.fit(dm, method='map', keys=keys, verbose=1)
        #model.fit(dm, method='norm_approx', keys=keys, verbose=1)
        model.fit(dm, method='mcmc', keys=keys, iter=100, thin=10, burn=1000, verbose=1)

        # remove all keys that are not relevant current model
        for k in dm.params.keys():
            if type(dm.params[k]) == dict:
                for j in dm.params[k].keys():
                    if not j in keys:
                        dm.params[k].pop(j)

    # post results to dismod_data_server
    url = dismod3.post_disease_model(dm)

    # form url to view results
    if opts.sex and opts.year and opts.region:
        url += '/%s/%s/%s' % (opts.region, opts.year, opts.sex)
    elif opts.region:
        url += '/%s' % opts.region

    # announce completion, and url to view results
    #tweet('%s fit complete %s' % (fit_str, url))
    sys.stdout.flush()

    # update job status file
    if opts.log and not opts.no_fit:
        if opts.type and not (opts.region and opts.sex and opts.year):
            dismod3.log_job_status(id, 'empirical_priors', opts.type, 'Completed')
        elif opts.region and opts.sex and opts.year and not opts.type:
            dismod3.log_job_status(id, 'posterior', '%s--%s--%s' % (opts.region, opts.sex, opts.year), 'Completed')
        elif opts.region and not (opts.sex and opts.year and opts.type):
            dismod3.log_job_status(id, 'within_each_region', opts.region, 'Completed')
        elif not (opts.region and opts.sex and opts.year and opts.type):
            dismod3.log_job_status(id, 'across_all_regions', 'all_regions', 'Completed')

if __name__ == '__main__':
    main()
