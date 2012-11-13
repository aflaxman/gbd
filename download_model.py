#!/usr/bin/python2.5
""" Copy model from web to j drive

Example
-------

$ ./run_on_cluster download_model.py 32142

"""

import optparse
import subprocess

import dismod3

def download_model(id):
    """ Copy model from web to j drive

    Parameters
    ----------
    id : int
      The model id number to copy
    """
    dir = dismod3.settings.JOB_WORKING_DIR % id  # TODO: refactor into a function

    try:
        model = dismod3.data.ModelData.load(dir)
        print 'model already on j drive in %s' % dir

    except (IOError, AssertionError):
        print 'downloading disease model'
        dm = dismod3.load_disease_model(id)

        import simplejson as json
        try:
            model = dismod3.data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))
        except Exception as e:
            print e
            print 'attempting to use old covariate format'
            import old_cov_data
            model = old_cov_data.from_gbd_jsons(json.loads(dm.to_json()))

        model.save(dir)
        print 'loaded data from json, saved in new format for next time in %s' % dir

if __name__ == '__main__':
    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    download_model(id)
