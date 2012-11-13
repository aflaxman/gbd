""" Remove all jobs from cluster

Example
-------

$ python kill_fit.py 4222

"""

import optparse
import subprocess

def kill_fit(id):
    call_str = 'qstat | grep %s | awk {\'print $1\'} | xargs qdel' % id
    print call_str
    subprocess.call(call_str, shell=True)

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

    kill_fit(id)


