#!/usr/local/bin/python2.5

"""
wrapper to run DisMod using command-line specified options
"""

import sys
import getopt

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def usage():
    print 'coming soon...'

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error, msg:
            raise Usage(msg)

        for o, a in opts:
            if o in ("-h", "--help"):
                usage()
                return 0
            
        try:
            dm_id = int(args[0])
        except IndexError:
            raise Usage('missing disease model id')

        #######################################################
        #######################################################
        # DisMod code goes here
        
        import dismod3
        import dismod3.generic_disease_model as probabilistic_model
        
        dm = dismod3.get_disease_model(dm_id)

        # merge in all-cause mortality rates from model 824
        # mort = dismod3.get_disease_model(824)
        # dm.data += mort.filter_data(data_type='all-cause mortality data')
        # dismod3.post_disease_model(dm)
        
        dm.set_param_age_mesh(range(30,91,2))
        dm.set_estimate_age_mesh(range(30,91))

        dm.set_priors('incidence data', 'smooth 10')
        dm.set_priors('prevalence data', 'smooth 10')
        
        probabilistic_model.initialize(dm)
        vars = probabilistic_model.setup(dm)
        map = probabilistic_model.map_fit(dm, vars)
        mc = probabilistic_model.mcmc_fit(dm, vars)
        url = dismod3.post_disease_model(dm)
        print url
        
        return 0

        #
        #######################################################
        #######################################################


    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2

if __name__ == "__main__":
    sys.exit(main())
