#!/usr/bin/python2.5
""" Distribute the parameter estimation throughout our cluster.

gbd_fit has options to start in daemon mode or in fitting mode.

Examples
--------

$ python computational_engine_daemon
"""

import optparse
import time
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

def log(message):
    sys.stdout.write('%s  gbd_fit  %s\n' % (time.strftime("%Y-%m-%d %H:%M:%S"), message))
    sys.stdout.flush()

def term(self, *args):
    log('dismod3 daemon received SIGTERM')
    sys.exit()

def main():
    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-d', '--daemon', default=True,
                      action='store_true', dest='daemon')

    parser.add_option('-l', '--log',
                      action='store_true', dest='log',
                      help='log the job running status')

    (options, args) = parser.parse_args()

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
            if not os.path.exists(dir):
                os.makedirs(dir)

            estimate_type = dm.params.get('run_status', {}).get('estimate_type', 'fit all individually')

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
                            time.sleep(1.)

            elif estimate_type.find('empirical priors') != -1:
                # fit empirical priors (by pooling data from all regions
                d = '%s/empirical_priors' % dir
                if os.path.exists(d):
                    rmtree(d)
                os.mkdir(d)
                os.mkdir('%s/stdout' % d)
                os.mkdir('%s/stderr' % d)
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

if __name__ == '__main__':
    main()
