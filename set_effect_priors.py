#!/usr/bin/python2.5
""" Set effect priors for model on j drive

Example
-------

$ ./run_on_cluster download_model.py 32142
$ ./run_on_cluster set_effect_priors.py 32142

"""

import optparse
import subprocess

import dismod3

def set_effect_priors(id):
    model_path = '/home/j/Project/dismod/output/dm-%d/'%id

    model = dismod3.data.load(model_path)

    for t in model.parameters:
        if 'fixed_effects' not in model.parameters[t]:
            continue

        for cov in model.input_data.filter(regex='^x_').columns:
            # set prior informative priors on all fixed effects,
            # to say that 10% effect is expected, and >30% is surprising
            model.parameters[t]['fixed_effects'][cov] = dict(dist='Normal', mu=0., sigma=.1)
            
    model.save(model_path)
    return model
    

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

    set_effect_priors(id)
