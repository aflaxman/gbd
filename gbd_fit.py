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

import dismod3
from dismod3.utils import clean, gbd_keys, type_region_year_sex_from_key
from dismod3.plotting import GBDDataHash

import sys
from os import popen

def tweet(message,
          user=dismod3.settings.DISMOD_TWITTER_NAME,
          password=dismod3.settings.DISMOD_TWITTER_PASSWORD):
    print 'tweeting %s for %s' % (message, user)

    message = '#dismod %s' % message
 
    url = 'http://twitter.com/statuses/update.xml' 
    curl = 'curl -s -u %s:%s -d status="%s" %s' % (user,password,message,url)
    try:
        pipe = popen(curl, 'r')
    except:
        pass

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

    (options, args) = parser.parse_args()

    if len(args) == 0:
        if options.daemon:
            try:
                tweet('starting dismod3 daemon...')
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
    while True:
        try:
            job_queue = dismod3.get_job_queue()
        except:
            job_queue = []
            
        for id in job_queue:
            tweet('processing job %d' % id)
            dismod3.remove_from_job_queue(id)
            dm = dismod3.get_disease_model(id)

            estimate_type = dm.params.get('estimate_type', 'fit all individually')

            # sort the regions so that the data rich regions are fit first
            data_hash = GBDDataHash(dm.data)
            sorted_regions = sorted(dismod3.gbd_regions, reverse=True,
                                    key=lambda r: len(data_hash.get(region=r)))
            
            if estimate_type.find('individually') != -1:
                #fit each region/year/sex individually for this model (84 processes!)
                for r in sorted_regions:
                    for s in dismod3.gbd_sexes:
                        for y in dismod3.gbd_years:
                            call_str = dismod3.settings.GBD_FIT_STR \
                                % ('-r %s -s %s -y %d' % (clean(r), s, y), id)
                            subprocess.call(call_str,
                                            shell=True)

            elif estimate_type.find('within each region') != -1:
                # fit each region individually, but borrow strength within gbd regions
                for r in sorted_regions:
                    subprocess.call(dismod3.settings.GBD_FIT_STR
                                    % ('-r %s' % clean(r), id), shell=True)

            elif estimate_type.find('across all regions') != -1:
                # fit all regions, years, and sexes together
                subprocess.call(dismod3.settings.GBD_FIT_STR % ('', id), shell=True)

            elif estimate_type.find('empirical priors') != -1:
                # fit empirical priors (by pooling data from all regions
                for t in ['case-fatality', 'remission', 'incidence', 'prevalence']:
                    subprocess.call(dismod3.settings.GBD_FIT_STR
                                    % ('-t %s' % t, id), shell=True)
                    
            else:
                tweet('unrecognized estimate type: %s' % estimate_type)

        time.sleep(dismod3.settings.SLEEP_SECS)
        
def fit(id, opts):
    fit_str = '(%d) %s %s %s' % (id, opts.region or '', opts.sex or '', opts.year or '')
    tweet('fitting disease model %s' % fit_str)

    dm = dismod3.get_disease_model(id)
    fit_str = '%s %s' % (dm.params['condition'], fit_str)

    sex_list = opts.sex and [ opts.sex ] or dismod3.gbd_sexes
    year_list = opts.year and [ opts.year ] or dismod3.gbd_years
    region_list = opts.region and [ opts.region ] or dismod3.gbd_regions
    keys = gbd_keys(region_list=region_list, year_list=year_list, sex_list=sex_list)
    if not opts.region:
        keys += gbd_keys(region_list=['world'], year_list=['total'], sex_list=['total'])

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
        fit_str += ' emp prior'
        import dismod3.logit_normal_model as model
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

        # fit the model with a normal approximation
        model.fit(dm, method='norm_approx', keys=keys, verbose=1)
        #model.fit(dm, method='map', keys=keys, verbose=1)

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
    tweet('%s fit complete %s' % (fit_str, url))
    
        
if __name__ == '__main__':
    main()
