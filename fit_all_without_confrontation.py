#!/usr/bin/python2.5
""" Fit all model parameters on cluster without data confrontation
"""

import optparse
import subprocess

import dismod3

def fit_all(id):
    """ Enqueues all jobs necessary to fit specified model
    to the cluster

    Parameters
    ----------
    id : int
      The model id number for the job to fit
    """

    if dismod3.settings.ON_SGE:
        print 'downloading disease model'
        dismod3.disease_json.create_disease_model_dir(id)
    dm = dismod3.load_disease_model(id)

    dir = dismod3.settings.JOB_WORKING_DIR % id  # TODO: refactor into a function

    # directory to save the country level posterior csv files
    temp_dir = dir + '/posterior/country_level_posterior_dm-' + str(id) + '/'

    #fit each region/year/sex individually for this model
    post_names = []
    for ii, r in enumerate(dismod3.gbd_regions):
        for s in dismod3.gbd_sexes:
            for y in dismod3.gbd_years:
                k = '%s+%s+%s' % (dismod3.utils.clean(r), dismod3.utils.clean(s), y)
                o = '%s/posterior/stdout/%s' % (dir, k)
                e = '%s/posterior/stderr/%s' % (dir, k)
                name_str = '%s%d%s%s%d' % (r[0], ii+1, s[0], str(y)[-1], id)
                post_names.append(name_str)

                if dismod3.settings.ON_SGE:
                    call_str = 'qsub -cwd -o %s -e %s ' % (o,e) \
                        + '-N %s ' % name_str \
                        + 'run_on_cluster.sh '
                else:
                    call_str = 'python '
                call_str += 'fit_without_confrontation.py %d -r %s -s %s -y %s' % (id, dismod3.utils.clean(r), dismod3.utils.clean(s), y)
                subprocess.call(call_str, shell=True)

    # after all posteriors have finished running, upload disease model json
    hold_str = '-hold_jid %s ' % ','.join(post_names)
    o = '%s/upload.stdout' % dir
    e = '%s/upload.stderr' % dir
    if dismod3.settings.ON_SGE:
        call_str = 'qsub -cwd -o %s -e %s ' % (o,e) \
            + hold_str \
            + '-N upld-%s ' % id \
            + 'run_on_cluster.sh '
    else:
        call_str = 'python '
    call_str += 'upload_fits.py %d' % id
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

    fit_all(id)


if __name__ == '__main__':
    main()
