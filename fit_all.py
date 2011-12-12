#!/usr/bin/python2.5
""" Fit all model parameters on cluster using empirical bayes.

Example
-------

$ python fit_all.py 4222    # submit jobs to cluster to estimate empirical priors followed by posteriors for model #4222

"""

import optparse
import subprocess

import dismod3
import data
reload(data)

def fit_all(id, consistent_empirical_prior=True, consistent_posterior=True, posteriors_only=False, posterior_types='pir', fast=False):
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
    dir = dismod3.settings.JOB_WORKING_DIR % id  # TODO: refactor into a function

    try:
        model = data.ModelData.load(dir)
        print 'loaded data from new format from %s' % dir

        # if we make it here, this model has already been run, so clean out the stdout/stderr dirs to make room for fresh messages
        call_str = 'rm -rf %s/empirical_priors/stdout/* %s/empirical_priors/stderr/* %s/posterior/stdout/* %s/posterior/stderr/* %s/json/dm-*-*.json' % (dir, dir, dir, dir, dir)
        print call_str
        subprocess.call(call_str, shell=True)

        # now load just the model, all previous fits are deleted
        dm = dismod3.load_disease_model(id)

    except (IOError, AssertionError):
        print 'downloading disease model'
        dm = dismod3.load_disease_model(id)

        import simplejson as json
        model = data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))
        model.save(dir)
        print 'loaded data from json, saved in new format for next time in %s' % dir

    o = '%s/empirical_priors/stdout/%d_running.txt' % (dir, id)
    f = open(o, 'w')
    import time
    f.write('Enqueued model %d on cluster at %s' % (id, time.strftime('%c')))
    f.close()



    # fit empirical priors (by pooling data from all regions)
    emp_names = []

    if not posteriors_only:

        if consistent_empirical_prior:
            t = 'all'
            o = '%s/empirical_priors/stdout/dismod_log_%s' % (dir, t)
            e = '%s/empirical_priors/stderr/dismod_log_%s' % (dir, t)
            name_str = '%s-%d' %(t[0], id)
            emp_names.append(name_str)
            if dismod3.settings.ON_SGE:
                call_str = 'qsub -cwd -o %s -e %s ' % (o, e) \
                    + '-N %s ' % name_str \
                    + 'run_on_cluster.sh '
            else:
                call_str = 'python '
            call_str += 'fit_world.py %d' % id

            if fast:
                call_str += ' --fast=true'
            
            subprocess.call(call_str, shell=True)

        else:
            for t in ['excess-mortality', 'remission', 'incidence', 'prevalence', 'prevalence_x_excess-mortality']:
                o = '%s/empirical_priors/stdout/dismod_log_%s' % (dir, t)
                e = '%s/empirical_priors/stderr/dismod_log_%s' % (dir, t)
                name_str = '%s-%d' %(t[0], id)
                emp_names.append(name_str)
                if dismod3.settings.ON_SGE:
                    call_str = 'qsub -cwd -o %s -e %s ' % (o, e) \
                        + '-N %s ' % name_str \
                        + 'run_on_cluster.sh '
                else:
                    call_str = 'python '
                call_str += 'fit_emp_prior.py %d -t %s' % (id, t)
                subprocess.call(call_str, shell=True)

    # directory to save the country level posterior csv files
    temp_dir = dir + '/posterior/country_level_posterior_dm-' + str(id) + '/'

    #fit each region/year/sex individually for this model
    hold_str = '-hold_jid %s ' % ','.join(emp_names)
    if posteriors_only:
        hold_str = ''
    post_names = []
    for ii, r in enumerate(dismod3.gbd_regions):
        for s in dismod3.gbd_sexes:
            for y in dismod3.gbd_years:
                k = '%s+%s+%s' % (dismod3.utils.clean(r), dismod3.utils.clean(s), y)
                o = '%s/posterior/stdout/dismod_log_%s' % (dir, k)
                e = '%s/posterior/stderr/dismod_log_%s' % (dir, k)
                name_str = '%s%d%s%s%d' % (r[0], ii+1, s[0], str(y)[-1], id)
                post_names.append(name_str)

                if dismod3.settings.ON_SGE:
                    call_str = 'qsub -cwd -o %s -e %s ' % (o,e) \
                        + hold_str \
                        + '-N %s ' % name_str \
                        + 'run_on_cluster.sh '
                else:
                    call_str = 'python '
                call_str += 'fit_posterior.py %d -r %s -s %s -y %s' % (id, dismod3.utils.clean(r), dismod3.utils.clean(s), y)

                if not consistent_posterior:
                    call_str += ' --inconsistent=True --types=%s' % posterior_types

                if fast:
                    call_str += ' --fast=true'

                subprocess.call(call_str, shell=True)

    # after all posteriors have finished running, upload disease model json
    hold_str = '-hold_jid %s ' % ','.join(post_names)
    o = '%s/empirical_priors/stdout/%d_upload.txt' % (dir, id)
    e = '%s/empirical_priors/stderr/%d_upload.txt' % (dir, id)
    if dismod3.settings.ON_SGE:
        call_str = 'qsub -cwd -o %s -e %s ' % (o,e) \
            + hold_str \
            + '-N upld-%s ' % id \
            + 'run_on_cluster.sh '
    else:
        call_str = 'python '
    call_str += 'upload_fits.py %d' % id
    subprocess.call(call_str, shell=True)

    return dm

def main():
    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-c', '--priorconsistent', default='True',
                      help='use consistent model for empirical priors')
    parser.add_option('-C', '--posteriorconsistent', default='True',
                      help='use consistent model for posteriors')
    parser.add_option('-t', '--posteriortypes', default='pir',
                      help='use consistent model for posteriors')
    parser.add_option('-o', '--onlyposterior', default='False',
                      help='skip empirical prior phase')
    parser.add_option('-f', '--fast', default='False',
                      help='use MAP only')
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    dm = fit_all(id,
                 consistent_empirical_prior=(options.priorconsistent.lower()=='true'),
                 consistent_posterior=(options.posteriorconsistent.lower()=='true'),
                 posteriors_only=(options.onlyposterior.lower()=='true'),
                 posterior_types=options.posteriortypes,
                 fast=(options.fast.lower() == 'true'))

if __name__ == '__main__':
    dm = main()
