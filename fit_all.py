""" Fit all model parameters on the cluster

Example
-------

$ python fit_all.py 14464    # submit jobs to cluster

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
    >>> fit_all.fit_all(14464)
    """

    # download the disease model json and store it in the working dir
    print 'downloading disease model'
    dismod3.disease_json.create_disease_model_dir(id)
    dm = dismod3.fetch_disease_model(id)
    
    # get the all-cause mortality data, and merge it into the model
    mort = dismod3.fetch_disease_model('all-cause_mortality')
    dm.data += mort.data
    dm.save()

    #fit each region/year/sex individually for this model
    post_names = []
    for ii, r in enumerate(dismod3.gbd_regions):
        for s in dismod3.gbd_sexes:
            for y in dismod3.gbd_years:
                dir = dismod3.settings.JOB_WORKING_DIR % id  # TODO: refactor into a function
                k = '%s.txt' % dismod3.utils.gbd_key_for('all', r, s, y)
                
                o = dir + '/posterior/stdout/%s' % k
                e = dir + '/posterior/stderr/%s' % k
                
                name_str = '%s%d%s%s%d' % (r[0], ii+1, s[0], str(y)[-1], id)
                post_names.append(name_str)


                if dismod3.settings.ON_SGE:
                    call_str = 'qsub -P ihme_dismod -cwd -o %s -e %s ' % (o,e) \
                           + '-N %s ' % name_str \
                           + 'run_on_cluster.sh '
                else:
                    call_str = 'python '
                call_str += 'fit_posterior.py %d -r %s -s %s -y %s' % (id, clean(r), s, y)
                subprocess.call(call_str, shell=True)

    # after all posteriors have finished running, upload disease model json
    hold_str = '-hold_jid %s ' % ','.join(post_names)
    o = dir + '/upload_stdout.txt'
    e = dir + '/upload_stderr.txt'

    if dismod3.settings.ON_SGE:
        call_str = 'qsub -P ihme_dismod -cwd -o %s -e %s ' % (o,e) \
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
