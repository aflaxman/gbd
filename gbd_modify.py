#!/usr/bin/python2.5
""" Command-line tool to change parameters for a DisMod model

Examples
--------

$ python2.5 gbd_modify.py 176 --ymax .3   # set the maximum y value on summary plots
$ python2.5 gbd_modify.py 176 --condition 'Type II Diabetes'  # set the condition that is being modelled
$ python2.5 gbd_modify.py 176 --notes 'testing covariate model' # set the notes for the model
$ python2.5 gbd_modify.py 176 --clone # create a clone of model 176

"""

import time
import optparse
import dismod3


def main():
    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)

    # flags for parameters to modify
    parser.add_option('-y', '--ymax',
                      help='set the maximum y value for summary plots')
    parser.add_option('-t', '--condition',
                      help='set the condition of the model')
    parser.add_option('-n', '--notes',
                      help='set the notes of the model')

    # boolean flags
    parser.add_option('-c', '--clone',
                      action='store_true', dest='clone',
                      help='create a clone of the model (leave specified model unchanged)')

    (opts, args) = parser.parse_args()

    # check that args are correct
    if len(args) == 1:
        try:
            id = int(args[0])
        except ValueError:
            parser.error('disease_model_id must be an integer')
            return
    else:
        parser.error('incorrect number of arguments')
        return

    # fetch requested model
    dm = dismod3.get_disease_model(id)

    # change values specified
    if opts.ymax:
        dm.set_ymax(float(opts.ymax))

    # TODO: get condition to actually change
    if opts.condition:
        dm.set_condition(opts.condition)
    if opts.notes:
        dm.set_notes(opts.notes)

    # clone if requested
    if opts.clone:
        dm.params.pop('id')  # dismod_data_server creates new model if it doesn't find an id

    # post results to dismod_data_server
    url = dismod3.post_disease_model(dm)

    # announce url to view results
    print url
    
        
if __name__ == '__main__':
    main()
