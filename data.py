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

        self.hierarchy = nx.DiGraph()
        self.hierarchy.add_node('all')

        self.nodes_to_fit = self.hierarchy.nodes()

    def save(self, path):
        """ Saves all model data in human-readable files

        Parameters
        ----------
        path : str, directory to save in

        Results
        -------
        Saves files to specified path, overwritting what was there
        before
        """

        pl.rec2csv(self.input_data.to_records(), path + 'input_data.csv')  # TODO: patch Pandas so that pandas.read_csv works when fields have commas in them
        pl.rec2csv(self.output_template.to_records(), path + 'output_template.csv')
        json.dump(self.parameters, open(path + 'parameters.json', 'w'), indent=2)
        json.dump(dict(nodes=sorted(self.hierarchy.nodes()), edges=sorted(self.hierarchy.edges())),
                  open(path + 'hierarchy.json', 'w'), indent=2)
        json.dump(self.nodes_to_fit, open(path + 'nodes_to_fit.json', 'w'), indent=2)

    @staticmethod
    def load(path):
        d = ModelData()

        d.input_data = pandas.DataFrame.from_records(pl.csv2rec(path + 'input_data.csv')).drop(['index'], 1) # TODO: patch Pandas so that pandas.read_csv works with pandas.DataFrame.to_csv
        d.output_template = pandas.DataFrame.from_records(pl.csv2rec(path + 'output_template.csv')).drop(['index'], 1)
        d.parameters = json.load(open(path + 'parameters.json'))

        hierarchy = json.load(open(path + 'hierarchy.json'))
        d.hierarchy.add_nodes_from(hierarchy['nodes'])
        d.hierarchy.add_edges_from(hierarchy['edges'])

        d.nodes_to_fit = json.load(open(path + 'nodes_to_fit.json'))

        return d

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
        dm = json.load(open(fname))

        # load some ancillary data from the gbd
        import dismod3
        import csv
        dm['countries_for'] = dict(
            [[dismod3.utils.clean(x[0]), x[1:]] for x in csv.reader(open(dismod3.settings.CSV_PATH + 'country_region.csv'))]
            )
        dm['population_by_age'] = dict(
            [[(r['Country Code'], r['Year'], r['Sex']),
              [max(.001,float(r['Age %d Population' % i])) for i in range(dismod3.settings.MAX_AGE)]] 
             for r in csv.DictReader(open(dismod3.settings.CSV_PATH + 'population.csv'))
             if len(r['Country Code']) == 3]
            )


        d = ModelData()

        d.input_data = ModelData._input_data_from_gbd_json(dm)
        d.output_template = ModelData._output_template_from_gbd_json(dm)
        d.parameters = ModelData._parameters_from_gbd_json(dm)
        d.hierarchy, d.nodes_to_fit = ModelData._hierarchy_from_gbd_json(dm)

        return d


    @staticmethod
    def _input_data_from_gbd_json(dm):
        """ translate input data"""
        import dismod3

        input_data = {}
        for field in 'data_type value sex age_start age_end year_start year_end standard_error effective_sample_size lower_ci upper_ci'.split():
            input_data[field] = [row.get(field) for row in dm['data']]
        input_data['area'] = [(row['country_iso3_code'] == 'all') and row['country_iso3_code'] or row['region'] for row in dm['data']]  # iso3 code or gbd region if iso3 code is blank or 'all'
        input_data['age_weights'] = [json.dumps(row['age_weights']) for row in dm['data']]  # store age_weights as json, since Pandas doesn't like arrays in arrays

        # add selected covariates
        for level in ['Country_level', 'Study_level']:
            for cv in dm['params']['covariates'][level]:
                if dm['params']['covariates'][level][cv]['rate']['value']:
                    input_data['x_%s'%cv] = [row.get(dismod3.utils.clean(cv)) for row in dm['data']]
        
        return pandas.DataFrame(input_data)


    @staticmethod
    def _output_template_from_gbd_json(dm):
        """ generate output template"""
        import dismod3
        output_template = {}
        for field in 'data_type area sex age_start age_end year_start year_end age_weights'.split():
            output_template[field] = []
        for level in ['Country_level', 'Study_level']:
            for cv in dm['params']['covariates'][level]:
                if dm['params']['covariates'][level][cv]['rate']['value']:
                    output_template['x_%s'%cv] = []

        for data_type in dismod3.settings.output_data_types:
            for region in dismod3.settings.gbd_regions[:3]:  # FIXME: just take first 3 regions, for faster development and testing
                for area in dm['countries_for'][dismod3.utils.clean(region)]:
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

                                age_weights = pl.array(dm['population_by_age'][area, year, sex][age_start:age_end])
                                age_weights = list(age_weights / age_weights.sum())
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
                                                
        return pandas.DataFrame(output_template)


    @staticmethod
    def _parameters_from_gbd_json(dm):
        """ copy expert priors"""
        parameters = ModelData().parameters
        old_name = dict(i='incidence', p='prevalence', rr='relative_risk', r='remission', f='excess_mortality', X='duration')
        for t in parameters:
            parameters[t]['parameter_age_mesh'] = dm['params']['global_priors']['parameter_age_mesh']
            parameters[t]['y_maximum'] = dm['params']['global_priors']['y_maximum']
            for prior in 'smoothness heterogeneity level_value level_bounds increasing decreasing'.split():
                parameters[t][prior] = dm['params']['global_priors'][prior][old_name[t]]
        return parameters


    @staticmethod
    def _hierarchy_from_gbd_json(dm):
        """ setup hierarchy and nodes_to_fit"""
        import dismod3

        superregions = [[15, 5, 9, 0, 12], [7, 8, 1], [17, 18, 19, 20], [14], [3], [4, 2, 16], [10, 11, 13, 6]]
        hierarchy = nx.DiGraph()
        nodes_to_fit = ['all']
        for i, superregion in enumerate(superregions):
            hierarchy.add_edge('all', 'super-region_%d'%i)
            for j in superregion:
                for year in [1990, 2005, 2010]:
                    for sex in 'male female'.split():
                        region = str(dismod3.utils.clean(dismod3.settings.gbd_regions[j]))
                        region_node = '%s+%d+%s' % (region, year, sex)
                        nodes_to_fit.append(region_node)
                        hierarchy.add_edge('super-region_%d'%i, region_node)
                        for iso3 in dm['countries_for'][region]:
                            hierarchy.add_edge(region_node, '%s+%d+%s'%(iso3,year,sex))

        return hierarchy, nodes_to_fit
