""" Data Handling Class for DisMod III"""

import pandas
import networkx as nx
import pylab as pl
import simplejson as json

class ModelData:
    """ ModelData object contains all information for a disease model:
        Data, model parameters, information about output
    """

    def __init__(self):
        self.input_data = pandas.DataFrame(columns=('data_type value area sex age_start age_end year_start year_end' +
                                           ' standard_error effective_sample_size lower_ci upper_ci age_weights').split())
        self.output_template = pandas.DataFrame(columns='data_type area sex age_start age_end year_start year_end age_weights'.split())
        self.parameters = dict(i={}, p={}, r={}, f={}, rr={}, X={})

        self.areas_hierarchy = nx.DiGraph()
        self.areas_hierarchy.add_node('all')

        self.areas_to_fit = self.areas_hierarchy.nodes()

        
    @staticmethod
    def from_gbd_json(fname):
        """ Create ModelData object from old DM3 JSON file

        Parameters
        ----------
        fname : str, filename of JSON file

        Results
        -------
        returns new ModelData object
        """
        import dismod3

        dm = json.load(open(fname))
        d = ModelData()

        # translate input data
        input_data = {}
        for field in 'data_type value sex age_start age_end year_start year_end standard_error effective_sample_size lower_ci upper_ci'.split():
            input_data[field] = [row.get(field) for row in dm['data']]
        input_data['area'] = [(row['region'] == 'all') and row['region'] or row['gbd_region'] for row in dm['data']]  # iso3 code or gbd region if iso3 code is blank or 'all'
        input_data['age_weights'] = [json.dumps(row['age_weights']) for row in dm['data']]  # store age_weights as json, since Pandas doesn't like arrays in arrays

        # add selected covariates
        for level in ['Country_level', 'Study_level']:
            for cv in dm['params']['covariates'][level]:
                if dm['params']['covariates'][level][cv]['rate']['value']:
                    input_data['x_%s'%cv] = [row.get(dismod3.utils.clean(cv)) for row in dm['data']]

        d.input_data = pandas.DataFrame(input_data)

        # generate output data
        import csv
        countries_for = dict(
            [[dismod3.utils.clean(x[0]), x[1:]] for x in csv.reader(open(dismod3.settings.CSV_PATH + 'country_region.csv'))]
            )
        output_template = {}
        for field in 'data_type area sex age_start age_end year_start year_end age_weights'.split():
            output_template[field] = []
        for level in ['Country_level', 'Study_level']:
            for cv in dm['params']['covariates'][level]:
                if dm['params']['covariates'][level][cv]['rate']['value']:
                    output_template['x_%s'%cv] = []

        for data_type in dismod3.settings.output_data_types:
            for region in dismod3.settings.gbd_regions[:3]:
                for area in countries_for[dismod3.utils.clean(region)]:
                    for year in dismod3.settings.gbd_years:
                        for sex in dismod3.settings.gbd_sexes:
                            sex = dismod3.utils.clean(sex)
                            for age_start, age_end in zip(dismod3.settings.gbd_ages[:-1], dismod3.settings.gbd_ages[1:]):
                                output_template['data_type'].append(data_type)
                                output_template['area'].append(area)
                                output_template['sex'].append(sex)
                                output_template['year_start'].append(float(year))
                                output_template['year_end'].append(float(year)+1)
                                output_template['age_start'].append(age_start)
                                output_template['age_end'].append(age_end)

                                age_weights = list(pl.ones(age_end-age_start)/float(age_end-age_start))  # TODO: get population age weights
                                output_template['age_weights'].append(json.dumps(list(age_weights)))

                                # merge in country level covariates
                                for level in ['Country_level', 'Study_level']:
                                    for cv in dm['params']['covariates'][level]:
                                        if dm['params']['covariates'][level][cv]['rate']['value']:
                                            if dm['params']['covariates'][level][cv]['value']['value'] == 'Country Specific Value':
                                                output_template['x_%s'%cv].append(dm['params']['derived_covariate'][cv].get('%s+%s+%s'%(area, year, sex)))
                                                if not output_template['x_%s'%cv]:
                                                    print 'WARNING: derived covariate %s not found for (%s, %s, %s)' % (cv, area, year, sex)

                                            else:
                                                output_template['x_%s'%cv].append(dm['params']['covariates'][level][cv]['value']['value'])
                                                
        d.output_template = pandas.DataFrame(output_template)

        # copy expert priors
        old_name = dict(i='incidence', p='prevalence', rr='relative_risk', r='remission', f='excess_mortality', X='duration')
        for t in d.parameters:
            d.parameters[t]['parameter_age_mesh'] = dm['params']['global_priors']['parameter_age_mesh']
            d.parameters[t]['y_maximum'] = dm['params']['global_priors']['y_maximum']
            for prior in 'smoothness heterogeneity level_value level_bounds increasing decreasing'.split():
                d.parameters[t][prior] = dm['params']['global_priors'][prior][old_name[t]]

        # setup areas hierarchy and areas_to_fit
        superregions = [[15, 5, 9, 0, 12], [7, 8, 1], [17, 18, 19, 20], [14], [3], [4, 2, 16], [10, 11, 13, 6]]
        d.areas_to_fit = ['all']
        for i, superregion in enumerate(superregions):
            d.areas_hierarchy.add_edge('all', 'super-region_%d'%i)
            for j in superregion:
                region = dismod3.utils.clean(dismod3.settings.gbd_regions[j])
                d.areas_to_fit.append(region)
                d.areas_hierarchy.add_edge('super-region_%d'%i, region)
                for iso3 in countries_for[region]:
                    d.areas_hierarchy.add_edge(region, iso3)

        return d
