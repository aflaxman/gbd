import pylab as pl
import pymc as mc

import dismod3

import generic_disease_model
import neg_binom_model

def setup(dm, keys):
    """ Generate the PyMC variables for a multi-region/year/sex generic
    disease model.

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
    
    Results
    -------
    vars : dict of PyMC stochs
      returns a dictionary of all the relevant PyMC objects for the
      multi-region/year/sex generic disease model.
    """
    
    vars = {}

    # for each region-year-sex triple among the keys
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                key = dismod3.utils.gbd_key_for('%s', r, y, s)
                if not key%'prevalence' in keys:
                    continue

                dm.set_units(key%'prevalence', '(per person)')
                dm.set_units(key%'duration', '(years)')
                for t in 'incidence', 'remission', 'excess-mortality':
                    dm.set_units(key%t, '(per person-year)')
                    #dm.get_initial_estimate(key%t, [d for d in dm.data if dm.relevant_to(d, t, r, y, s)])

                data = [d for d in dm.data if dm.relevant_to(d, 'all', r, y, s)]
                # Also include data that is not sex specific (i.e. sex == 'total') here
                data += [d for d in dm.data if dm.relevant_to(d, 'all', r, y, 'total')]  # FIXME: this will double include data with sex == 'all'
                
                sub_vars = generic_disease_model.setup(dm, key, data)
                vars.update(sub_vars)
    
    return vars
