#!/usr/bin/python2.5
""" Distribute the parameter estimation throughout our cluster.
"""

import time

import fit_all
import dismod3

def daemon_loop():
    while True:
        try:
            job_queue = dismod3.disease_json.get_job_queue()
        except:
            job_queue = []
        
        for param_id in job_queue:
            job_params = dismod3.disease_json.remove_from_job_queue(param_id)
            id = int(job_params['dm_id'])
            fit_all.fit_all(id)
        time.sleep(dismod3.settings.SLEEP_SECS)

if __name__ == '__main__':
    daemon_loop()
