import glob
import os
import shutil
import optparse

import dismod3

def move_rerun_to_new(old_id, new_id):
    dismod3.disease_json.JOB_WORKING_DIR = '/home/j/Project/dismod/dismod_status/prod/dm_best_rerun_2012_06_27/dm-%d'
    dm = dismod3.disease_json.load_disease_model(old_id)

    dm.id = new_id
    dm.params['id'] = new_id

    dismod3.disease_json.post_disease_model(dm)

def rename_posterior_draws_in_new(old_id, new_id):
    to_rename = glob.glob('/home/j/Project/dismod/dismod_status/prod/dm-%d/posterior/dm-%d*' % (new_id, old_id))
    #print to_rename
    for f in to_rename:
        #print f, 'becomes', f.replace(str(old_id), str(new_id))
        os.rename(f, f.replace(str(old_id), str(new_id)))
                                                          
if __name__ == '__main__':
    usage = 'usage: %prog [options] old_id new_id'
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.error('incorrect number of arguments')

    try:
        old_id = int(args[0])
        new_id = int(args[1])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    assert not os.path.exists('/home/j/Project/dismod/dismod_status/prod/dm-%d' % new_id)
    move_rerun_to_new(old_id, new_id)
    shutil.copytree('/home/j/Project/dismod/dismod_status/prod/dm_best_rerun_2012_06_27/dm-%d'%old_id,
                    '/home/j/Project/dismod/dismod_status/prod/dm-%d'%new_id)
    rename_posterior_draws_in_new(old_id, new_id)
