#!/usr/bin/python2.5
""" Fit all model parameters on cluster using empirical bayes.

Example
-------

$ ./fit_all 4222    # submit jobs to cluster to estimate empirical priors followed by posteriors for model #4222

"""

import optparse
import os
import subprocess
from shutil import rmtree

import dismod3
from dismod3.utils import clean, gbd_keys, type_region_year_sex_from_key

def main():
    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    if not dismod3.settings.ON_SGE:
        parser.error('dismod3.settings.ON_SGE must be true to fit_all automatically')

    #dm = dismod3.get_disease_model(id)

    # make directory structure to store computation output
    dir = dismod3.settings.JOB_WORKING_DIR % id
    if os.path.exists(dir):
        rmtree(dir)
    os.makedirs(dir)

    for phase in ['empirical_priors', 'posterior']:
        os.mkdir('%s/%s' % (dir, phase))
        for f_type in ['stdout', 'stderr']:
            os.mkdir('%s/%s/%s' % (dir, phase, f_type))


    # fit empirical priors (by pooling data from all regions
    emp_names = []
    for t in ['excess-mortality', 'remission', 'incidence', 'prevalence']:
        o = '%s/empirical_priors/stdout/%s' % (dir, t)
        e = '%s/empirical_priors/stderr/%s' % (dir, t)
        GBD_FIT_STR = 'qsub -cwd -o %s -e %s /home/OUTPOST/abie/gbd/gbd_fit.sh %s %d'
        name_str = '%s-%d' %(t[0], id)
        emp_names.append(name_str)
        call_str = 'qsub -cwd -o %s -e %s ' % (o, e) \
                        + '-N %s ' % name_str \
                        + 'run_on_cluster.sh /home/OUTPOST/abie/gbd_dev/gbd/fit_emp_prior.py %d -t %s' % (id, t)
        subprocess.call(call_str, shell=True)
    hold_str = '-hold_jid %s ' % ','.join(emp_names)
    #fit each region/year/sex individually for this model
    for ii, r in enumerate(dismod3.gbd_regions):
        for s in dismod3.gbd_sexes:
            for y in dismod3.gbd_years:
                k = '%s+%s+%s' % (clean(r), s, y)
                o = '%s/posterior/stdout/%s' % (dir, k)
                e = '%s/posterior/stderr/%s' % (dir, k)
                call_str = 'qsub -cwd -o %s -e %s ' % (o,e) \
                           + hold_str \
                           + '-N %s%d%s%s-%d ' % (r[0], ii+1, s[0], str(y)[-2:], id) \
                           + 'run_on_cluster.sh fit_posteriors.py %d -r %s -s %s -y %s' % (id, clean(r), s, y)
                subprocess.call(call_str, shell=True)
        
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
        if opts.type == 'all':
            for t in ['incidence', 'prevalence', 'remission', 'excess-mortality']:
                model.fit_emp_prior(dm, t)
        else:
            model.fit_emp_prior(dm, opts.type)

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
        model.fit(dm, method='map', keys=keys, verbose=1)
        model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1)
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
    if opts.sex and opts.year and opts.region:
        url += '/%s/%s/%s' % (opts.region, opts.year, opts.sex)
    elif opts.region:
        url += '/%s' % opts.region

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
