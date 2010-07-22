#!/usr/bin/python2.5
""" Generate empirical prior of specified parameter type

Expects the disase model json to be saved already.
"""

import simplejson as json
import dismod3

def upload_fits(id):
    """ Send results of cluster fits to dismod server

    Parameters
    ----------
    id : int
      The model id number

    Example
    -------
    >>> import fit_emp_prior
    >>> fit_emp_prior.fit_emp_prior(2552, 'incidence')
    >>> import upload_fits
    >>> upload_fits.upload_fits(2552)
    """
    # load disease model
    dm = dismod3.load_disease_model(id)  # this merges together results from all fits
    dismod3.try_posting_disease_model(dm, ntries=5)

def main():
    import optparse

    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    upload_fits(id)
      

if __name__ == '__main__':
    main()
