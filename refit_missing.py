#!/usr/bin/python2.5
""" Fit all model parameters on cluster using empirical bayes.

Example
-------

$ python fit_all.py 4222    # submit jobs to cluster to estimate empirical priors followed by posteriors for model #4222

"""

import optparse
import subprocess

import dismod3

def refit_missing(id, consistent_posterior=True, posterior_types='p i r', fast=False):
    """ Enqueues all jobs necessary to fit missing regions for specified model
    to the cluster

    Parameters
    ----------
    id : int
      The model id number for the job to fit
    """
    dir = dismod3.settings.JOB_WORKING_DIR % id  # TODO: refactor into a function

    import data
    model = data.ModelData.load(dir)

    o = '%s/empirical_priors/stdout/%d_running.txt' % (dir, id)
    f = open(o, 'a')
    import time
    f.write('\nEnqueued model %d on cluster at %s' % (id, time.strftime('%c')))
    f.close()

    # directory to save the country level posterior csv files
    temp_dir = dir + '/posterior/country_level_posterior_dm-' + str(id) + '/'

    #fit each region/year/sex individually for this model
    hold_str = ''
    post_names = []
    pretty_names = ''
    for ii, r in enumerate(dismod3.gbd_regions):
        for s in dismod3.gbd_sexes:
            for y in dismod3.gbd_years:
                k = '%s+%s+%s' % (dismod3.utils.clean(r), dismod3.utils.clean(s), y)
                o = '%s/posterior/stdout/dismod_log_%s' % (dir, k)
                e = '%s/posterior/stderr/dismod_log_%s' % (dir, k)
                name_str = '%s%d%s%s%d' % (r[0], ii+1, s[0], str(y)[-1], id)

                # if json file exists for this model, then continue
                try:
                    f = open('%s/json/dm-%d-posterior-%s-%s-%s.json'%(dir, id, dismod3.utils.clean(r), dismod3.utils.clean(s), y))
                    f.close()
                except IOError:
                    post_names.append(name_str)
                    pretty_names += 'http://winthrop.ihme.washington.edu/dismod/show/tile_%d_xxx+all+%s+%s+%s.png\n' % (id, dismod3.utils.clean(r), y, dismod3.utils.clean(s))

                    if dismod3.settings.ON_SGE:
                        #call_str = 'qsub -cwd -o %s -e %s ' % (o,e) \
                        call_str = 'qsub -cwd ' \
                            + hold_str \
                            + '-N %s ' % name_str \
                            + 'run_on_cluster.sh '
                    else:
                        call_str = 'python '
                    call_str += 'fit_posterior.py %d -r %s -s %s -y %s' % (id, dismod3.utils.clean(r), dismod3.utils.clean(s), y)

                    if not consistent_posterior:
                        call_str += ' --inconsistent=True --types="%s"' % posterior_types

                    if fast:
                        call_str += ' --fast=true'

                    subprocess.call(call_str, shell=True)

    # after all posteriors have finished running, upload disease model json
    if len(post_names) > 0:
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
    else:
        print 'Nothing found missing to refit'
        
    print pretty_names


def main():
    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-C', '--posteriorconsistent', default='True',
                      help='use consistent model for posteriors')
    parser.add_option('-t', '--posteriortypes', default='p i r',
                      help='use consistent model for posteriors')
    parser.add_option('-f', '--fast', default='False',
                      help='use MAP only')
    parser.add_option('-r', '--report', default='False',
                      help='report only')
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    if options.report.lower() == 'false':
        dm = refit_missing(id,
                           consistent_posterior=(options.posteriorconsistent.lower()=='true'),
                           posterior_types=options.posteriortypes,
                           fast=(options.fast.lower() == 'true'))
    else:
        report(id)
        

def report(id):
    """ Report all errors from stderr for id, and any missing files

    Parameters
    ----------
    id : int
      The model id number to report on
    """
    dir = dismod3.settings.JOB_WORKING_DIR % id  # TODO: refactor into a function

    import data
    model = data.ModelData.load(dir)

    o = '%s/empirical_priors/stdout/%d_running.txt' % (dir, id)
    f = open(o)
    print f.read()
    f.close()

    #fit each region/year/sex individually for this model
    for ii, r in enumerate(dismod3.gbd_regions):
        for s in dismod3.gbd_sexes:
            for y in dismod3.gbd_years:
                k = '%s+%s+%s' % (dismod3.utils.clean(r), dismod3.utils.clean(s), y)
                o = '%s/posterior/stdout/dismod_log_%s' % (dir, k)
                e = '%s/posterior/stderr/dismod_log_%s' % (dir, k)

                # if json file exists for this model, then continue
                try:
                    f = open('%s/json/dm-%d-posterior-%s-%s-%s.json'%(dir, id, dismod3.utils.clean(r), dismod3.utils.clean(s), y))
                    f.close()
                except IOError:
                    print '\n\nJSON not found for fit_posterior.py %d -r %s -s %s -y %s' % (id, dismod3.utils.clean(r), dismod3.utils.clean(s), y)
                    f = open(e)
                    print f.read()
                    f.close()
                    print
                    
if __name__ == '__main__':
    dm = main()
