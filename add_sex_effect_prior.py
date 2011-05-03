#!/usr/bin/python2.5
""" Add an "effect prior" to a given disease model

Example
-------

$ python add_sex_effect_prior.py 11094 prevalence 2.56 1.46 4.49  # set sex effect to this (mean, lower ci, upper ci) 

"""

import sys
import simplejson as json

import dismod3




def main():
    argv = sys.argv
    assert len(argv) > 2, 'missing options'

    id = int(argv[1])
    
    assert len(argv) == 6, 'usage: add_sex_effect_prior.py model_id rate_type mean lower_ci upper_ci'
    add_sex_effect_prior(id=int(argv[1]), type=argv[2],
                         mean=float(argv[3]),
                         lower_ci=float(argv[4]), upper_ci=float(argv[5]))

def add_sex_effect_prior(id, type, mean, lower_ci, upper_ci):
    dm = dismod3.disease_json.DiseaseJson(json.dumps({'params': {}, 'data': [], 'id': id}))
    dm.params['sex_effect_%s'%type] = dict(mean=mean, upper_ci=upper_ci, lower_ci=lower_ci)
    dismod3.post_disease_model(dm)

if __name__ == '__main__':
    main()
