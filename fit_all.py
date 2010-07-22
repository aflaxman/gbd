#!/usr/bin/python2.5
""" Fit all model parameters on cluster using empirical bayes.

Example
-------

$ python fit_all.py 4222    # submit jobs to cluster to estimate empirical priors followed by posteriors for model #4222

"""

import optparse
import os
import subprocess
from shutil import rmtree

import dismod3
from dismod3.utils import clean, gbd_keys, type_region_year_sex_from_key


def fit_all(id):
    """ Enqueues all jobs necessary to fit specified model
    to the cluster

    Parameters
    ----------
    id : int
      The model id number for the job to fit

    Example
    -------
    >>> import fit_all
    >>> fit_all.fit_all(2552)
    """
    # make directory structure to store computation output
    dir = dismod3.settings.JOB_WORKING_DIR % id
    if os.path.exists(dir):
        rmtree(dir)  # TODO: this fails when files are "busy", could move to dir + '_' or something instead
    os.makedirs(dir)

    for phase in ['empirical_priors', 'posterior']:
        os.mkdir('%s/%s' % (dir, phase))
        for f_type in ['stdout', 'stderr', 'pickle']:
            os.mkdir('%s/%s/%s' % (dir, phase, f_type))
    os.mkdir('%s/json' % dir)
    os.mkdir('%s/png' % dir)

    # download the disease model json and store it in the working dir
    print 'downloading disease model'
    dm = dismod3.fetch_disease_model(id)
    
    # get the all-cause mortality data, and merge it into the model
    mort = dismod3.fetch_disease_model('all-cause_mortality')
    dm.data += mort.data
    dm.save()

    # fit empirical priors (by pooling data from all regions)
    emp_names = []
    for t in ['excess-mortality', 'remission', 'incidence', 'prevalence']:
        o = '%s/empirical_priors/stdout/%s' % (dir, t)
        e = '%s/empirical_priors/stderr/%s' % (dir, t)
        name_str = '%s-%d' %(t[0], id)
        emp_names.append(name_str)
        call_str = 'qsub -cwd -o %s -e %s ' % (o, e) \
                        + '-N %s ' % name_str \
                        + 'run_on_cluster.sh /home/OUTPOST/abie/gbd_dev/gbd/fit_emp_prior.py %d -t %s' % (id, t)
        subprocess.call(call_str, shell=True)

    #fit each region/year/sex individually for this model
    hold_str = '-hold_jid %s ' % ','.join(emp_names)
    post_names = []
    for ii, r in enumerate(dismod3.gbd_regions):
        for s in dismod3.gbd_sexes:
            for y in dismod3.gbd_years:
                k = '%s+%s+%s' % (clean(r), s, y)
                o = '%s/posterior/stdout/%s' % (dir, k)
                e = '%s/posterior/stderr/%s' % (dir, k)
                name_str = '%s%d%s%s%d' % (r[0], ii+1, s[0], str(y)[-1], id)
                post_names.append(name_str)
                call_str = 'qsub -cwd -o %s -e %s ' % (o,e) \
                           + hold_str \
                           + '-N %s ' % name_str \
                           + 'run_on_cluster.sh fit_posterior.py %d -r %s -s %s -y %s' % (id, clean(r), s, y)
                subprocess.call(call_str, shell=True)

    # after all posteriors have finished running, upload disease model json
    hold_str = '-hold_jid %s ' % ','.join(post_names)
    o = '%s/upload.stdout' % dir
    e = '%s/upload.stderr' % dir
    call_str = 'qsub -cwd -o %s -e %s ' % (o,e) \
               + hold_str \
               + '-N upld-%s ' % id \
               + 'run_on_cluster.sh upload_fits.py %d' % id
    subprocess.call(call_str, shell=True)


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

    fit_all(id)


if __name__ == '__main__':
    main()
