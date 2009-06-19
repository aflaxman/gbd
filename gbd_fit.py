#!/usr/bin/python2.5
""" Distribute the parameter estimation throughout our cluster.

gbd_fit has options to start in daemon mode or in fitting mode.

Examples
--------

$ python gbd_fit --daemon    # launch daemon that will fit models as they become available
$ python gbd_fit 10   # launch fitting calculation to estimate parameters for model #10

"""

import time
import optparse
import subprocess

import dismod3
from dismod3.utils import clean, type_region_year_sex_from_key
from dismod3.plotting import GBDDataHash

import sys
from os import popen

def tweet(message,
          user=dismod3.settings.DISMOD_TWITTER_NAME,
          password=dismod3.settings.DISMOD_TWITTER_PASSWORD):
    print 'tweeting %s for %s' % (message, user)

    message += ' #dismod'
 
    url = 'http://twitter.com/statuses/update.xml' 
    curl = 'curl -s -u %s:%s -d status="%s" %s' % (user,password,message,url)

    pipe = popen(curl, 'r')

def main():
    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-s', '--sex', dest='sex',
                      help='only estimate given sex (valid settings ``male``, ``female``, ``all``)')
    parser.add_option('-y', '--year',
                      help='only estimate given year (valid settings ``1990``, ``2005``)')
    parser.add_option('-r', '--region',
                      help='only estimate given GBD Region')

    parser.add_option('-P', '--prevprior',
                      help='append string to prevalence priors')
    parser.add_option('-I', '--inciprior',
                      help='append string to incidence priors')
    parser.add_option('-C', '--caseprior',
                      help='append string to case-fatality priors')
    parser.add_option('-R', '--remiprior',
                      help='append string to remission priors')
    parser.add_option('-N', '--newprior',
                      action='store_true', dest='new_prior',
                      help='replace (instead of appending) prior strs')

    parser.add_option('-d', '--daemon',
                      action='store_true', dest='daemon')

    (options, args) = parser.parse_args()

    if len(args) == 0:
        if options.daemon:
            try:
                daemon_loop()
            finally:
                tweet('...dismod3 daemon shutting down')
        else:
            parser.error('incorrect number of arguments')
    elif len(args) == 1:
        try:
            id = int(args[0])
        except ValueError:
            parser.error('disease_model_id must be an integer')
        fit(id, options)
    else:
        parser.error('incorrect number of arguments')

def daemon_loop():
    tweet('starting dismod3 daemon...')

    while True:
        try:
            job_queue = dismod3.get_job_queue()
        except:
            job_queue = []
            
        for id in job_queue:
            tweet('processing job %d' % id)
            dismod3.remove_from_job_queue(id)
            dm = dismod3.get_disease_model(id)

            estimation_type = dm.params.get('estimation_type', 'fit all individually')
            
            if estimation_type.find('individually') != -1:
                #fit each region/year/sex individually for this model (84 processes!)
                data_hash = GBDDataHash(dm.data)
                sorted_regions = sorted(dismod3.gbd_regions, reverse=True,
                                        key=lambda r: len(data_hash.get(region=r)))
                for r in sorted_regions:
                    for s in dismod3.gbd_sexes:
                        for y in dismod3.gbd_years:
                            call_str = dismod3.settings.GBD_FIT_STR \
                                % ('-r %s -s %s -y %d' % (clean(r), s, y), id)
                            subprocess.call(call_str,
                                            shell=True)
            elif estimation_type.find('within regions') != -1:
                # fit each region individually, but borrow strength within gbd regions
                for r in dismod3.gbd_regions:
                    subprocess.call(dismod3.settings.GBD_FIT_STR
                                    % ('-r %s' % clean(r), id), shell=True)
            elif estimation_type.find('between regions') != -1:
                # fit all regions, years, and sexes together
                subprocess.call(dismod3.settings.GBD_FIT_STR % ('', id), shell=True)
            else:
                tweet('unrecognized estimation type: %s' % etimation_type)
        time.sleep(dismod3.settings.SLEEP_SECS)
        
def fit(id, opts):
    import dismod3.gbd_disease_model as model
    
    fit_str = '%d %s %s %s' % (id, opts.region, opts.sex, opts.year)
    tweet('fitting disease model %s' % fit_str)

    dm = dismod3.get_disease_model(id)

    # get the all-cause mortality data, and merge it into the model
    mort = dismod3.get_disease_model('all-cause_mortality')
    dm.data += mort.data

    sex_list = opts.sex and [ opts.sex ] or dismod3.gbd_sexes
    year_list = opts.year and [ opts.year ] or dismod3.gbd_years
    region_list = opts.region and [ opts.region ] or dismod3.gbd_regions

    # fit individually, if sex, year, and region are specified
    if opts.sex and opts.year and opts.region:
        dm.params['estimate_type'] = 'fit individually'
        
    keys = model.gbd_keys(region_list=region_list, year_list=year_list, sex_list=sex_list)
    
    if not opts.region:
        keys += model.gbd_keys(region_list=['world'], year_list=['total'], sex_list=['total'])

    # quick way to add/replace priors from the command line
    for rate_type, additional_priors in [ ['prevalence', opts.prevprior], ['incidence', opts.inciprior],
                                          ['remission', opts.remiprior], ['case-fatality', opts.caseprior] ]:
        if additional_priors:
            for k in keys:
                if  type_region_year_sex_from_key(k)[0] == rate_type:
                    if opts.new_prior:
                        dm.set_priors(k, additional_priors)
                    else:
                        dm.set_priors(k, dm.get_priors(k) + additional_priors)
    # TODO:  make sure that the post_disease_model only stores the parts we want it to
    model.fit(dm, method='map', keys=keys)

    # remove all keys that are not relevant current model
    for k in dm.params.keys():
        if type(dm.params[k]) == dict:
            for j in dm.params[k].keys():
                if not j in keys:
                    dm.params[k].pop(j)
                    
    url = dismod3.post_disease_model(dm)

    if opts.sex and opts.year and opts.region:
        url += '?sex=%s&year=%s&region=%s' % (opts.region, opts.year, opts.sex))

    tweet('initial fit of %s complete %s' % (fit_str, url))
                    
    model.fit(dm, method='mcmc', keys=keys)
    dismod3.post_disease_model(dm)
    tweet('MCMC fit of %s complete %s' % (fit_str, url))

    model.fit(dm, method='map', keys=keys)
    dismod3.post_disease_model(dm)
    tweet('final fit of %s' % (fit_str, url))
    
        
if __name__ == '__main__':
    main()
