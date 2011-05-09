#!/usr/bin/python2.5
""" Add an "effect prior" to a given disease model

Example
-------

$ python add_region_effect_prior.py 11094 uninformative

"""

import sys
import simplejson as json

import dismod3




def main():
    argv = sys.argv
    assert len(argv) == 3, 'missing options'

    id = int(argv[1])
    
    add_region_effect_prior(id=int(argv[1]), val=argv[2])

def add_region_effect_prior(id, val):
    dm = dismod3.disease_json.DiseaseJson(json.dumps({'params': {}, 'data': [], 'id': id}))
    dm.params['region_effects'] = val
    dismod3.post_disease_model(dm)

if __name__ == '__main__':
    main()
