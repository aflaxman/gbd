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
from dismod3.utils import clean

def main():
    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-s', '--sex', dest='sex',
                      help='only estimate given sex (valid settings ``male``, ``female``, ``all``)')
    parser.add_option('-y', '--year',
                      help='only estimate given year (valid settings ``1990``, ``2005``)')
    parser.add_option('-r', '--region',
                      help='only estimate given GBD Region')
    parser.add_option('-d', '--daemon',
                      action='store_true', dest='daemon')

    (options, args) = parser.parse_args()

    if len(args) == 0:
        if options.daemon:
            daemon_loop()
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
    print 'starting dismod3 daemon...'

    while True:
        job_queue = dismod3.get_job_queue()
        for id in job_queue:
            print 'processing job %d' % id
            dismod3.remove_from_job_queue(id)
            dm = dismod3.get_disease_model(id)

            estimation_type = dm.params.get('estimation_type', 'fit all individually')
            import pdb; pdb.set_trace()
            
            if estimation_type.find('individually') != -1:
                #fit each region/year/sex individually for this model (84 processes!)
                for r in dismod3.gbd_regions:
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
                                    % ('-r %s' % r, id))
            elif estimateion_type.find('between regions') != -1:
                # fit all regions, years, and sexes together
                subprocess.call(dismod3.settings.GBD_FIT_STR % ('', id))
            else:
                print 'unrecognized estimation type: %s' % etimation_type
        time.sleep(dismod3.settings.SLEEP_SECS)
        
def fit(id, opts):
    import dismod3.gbd_disease_model as model

    print 'fitting disease model %d' % id

    dm = dismod3.get_disease_model(id)

    # get the all-cause mortality data, and merge it into the model
    mort = dismod3.get_disease_model('all-cause_mortality')
    dm.data += mort.data

    #import pdb; pdb.set_trace()
    
    sex_list = opts.sex and [ opts.sex ] or dismod3.gbd_sexes
    year_list = opts.year and [ opts.year ] or dismod3.gbd_years
    region_list = opts.region and [ opts.region ] or dismod3.gbd_regions

    keys = model.gbd_keys(region_list=region_list, year_list=year_list, sex_list=sex_list)

    #if not opts.sex or not opts.year:
    #    for r in region_list:
    #        keys.append('%s-sex-year-similarity-potential' % r)
    
    if not opts.region:
        keys += model.gbd_keys(region_list=['world'], year_list=['total'], sex_list=['total'])
        #keys.append('world-similarity-potential')
    
    # TODO:  make sure that the post_disease_model only stores the parts we want it to
    model.fit(dm, method='map', keys=keys)

    # remove all keys that are not relevant current model
    for k in dm.params:
        if type(dm.params[k]) == dict:
            for j in dm.params[k]:
                if not j in keys:
                    dm.params[k].pop(j)
            
    dismod3.post_disease_model(dm)

    model.fit(dm, method='mcmc', keys=keys)
    dismod3.post_disease_model(dm)

    model.fit(dm, method='map', keys=keys)
    dismod3.post_disease_model(dm)
    
        
if __name__ == '__main__':
    main()
