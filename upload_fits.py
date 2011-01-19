#!/usr/bin/python2.5
""" Generate empirical prior of specified parameter type

Expects the disase model json to be saved already.
"""

import simplejson as json
import dismod3
import zipfile, os
from shutil import rmtree

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
    #zip_country_level_posterior_files(id)

def zip_country_level_posterior_files(id):
    """  Zip country level posterior files in the directory of the
    job_working_directory/posterior/country_level_posterior_dm-'id', 
    and then remove the directory containing the files

    Parameters
    ----------
    id : int
      The model id number
    """
    # job working directory
    job_wd = dismod3.settings.JOB_WORKING_DIR % id

    # directory containing the csv files
    directory = 'country_level_posterior_dm-' + str(id)

    try:
        # move to directory
        orig_dir = os.getcwd()
        os.chdir(job_wd + '/posterior/')

        # open an archive for writing
        a = zipfile.ZipFile(directory + '.zip', 'w', zipfile.ZIP_DEFLATED)

        # put files into the archive
        for f in os.listdir(directory):
            print "archiving file %s" % f
            a.write(directory + '/' + f)

        # close the archive
        a.close()

        # remove directory
        rmtree(directory)

        # move back directory
        os.chdir(orig_dir)

    except Exception,e:
        print e

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
