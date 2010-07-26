#!/usr/bin/python2.5
""" Poll the dismod data server for jobs to run, and call fit_all on them

Examples
--------

$ python fit_daemon

"""
import time

import dismod3
import fit_all

def daemon_loop():
    while True:
        try:
            job_queue = dismod3.get_job_queue()
        except:
            job_queue = []
        
        for param_id in job_queue:
            print 'processing job %d' % param_id
            job_params = dismod3.remove_from_job_queue(param_id)
            id = int(job_params['dm_id'])
            fit_all.fit_all(id)
            
        time.sleep(dismod3.settings.SLEEP_SECS)

def main():
    try:
        print 'starting dismod3 daemon...'
        daemon_loop()
    finally:
        print 'dismod3 daemon shutting down'

if __name__ == '__main__':
    main()
