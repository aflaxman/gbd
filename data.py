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
        self.output_template = pandas.DataFrame(columns='data_type area sex year pop'.split())
        self.parameters = dict(i={}, p={}, r={}, f={}, rr={}, X={}, ages=range(101))

        self.hierarchy = nx.DiGraph()
        self.hierarchy.add_node('all')

        self.nodes_to_fit = self.hierarchy.nodes()

    def get_data(self, data_type):
        if len(self.input_data) > 0:
            return self.input_data[self.input_data['data_type'] == data_type]
        else:
            return self.input_data


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
        json.dump(dict(nodes=[[n, self.hierarchy.node[n]] for n in sorted(self.hierarchy.nodes())],
                       edges=[[u, v, self.hierarchy.edge[u][v]] for u,v in sorted(self.hierarchy.edges())]),
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

        print 'loading %s' % fname
        dm = json.load(open(fname))
        return ModelData.from_gbd_jsons(dm)


    @staticmethod
    def from_gbd_jsons(dm):
        """ Create ModelData object from old DM3 JSON file

        Parameters
        ----------
        dm : str, the JSON data

        Results
        -------
        returns new ModelData object
        """
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

        print 'load completed successfully'

        return d


    @staticmethod
    def _input_data_from_gbd_json(dm):
        """ translate input data"""
        import dismod3

        input_data = {}
        for field in 'effective_sample_size age_start age_end year_start year_end'.split():
            input_data[field] = []
            for row in dm['data']:
                val = row.get(field, '')
                if val == '':
                    val = pl.nan
                input_data[field].append(float(val))

        input_data['sex'] = []
        for row in dm['data']:
            input_data['sex'].append(row['sex'])

            # replace sex 'all' with sex 'total'
            if input_data['sex'][-1] == 'all':
                input_data['sex'][-1] = 'total'
                
            assert input_data['sex'][-1] != ''

        new_type_name = {'incidence data':'i', 'prevalence data': 'p', 'remission data': 'r', 'excess-mortality data': 'f',
                         'prevalence x excess-mortality data': 'pf', 'all-cause mortality data': 'm', 'relative-risk data': 'rr',
                         'duration data': 'X', 'smr data': 'smr', 'cause-specific mortality data': 'csmr', 'mortality data': 'm_with'}
        input_data['data_type'] = [new_type_name[row['data_type']] for row in dm['data']]

        for field in 'value standard_error lower_ci upper_ci'.split():
            input_data[field] = [float(row.get(field, '') or pl.nan) / float(row.get('units', '1').replace(',', '')) for row in dm['data']]

        input_data['area'] = []
        for row in dm['data']:
            val = row.get('country_iso3_code', '')
            if val == '' or val == 'all':
                val = dismod3.utils.clean(row['gbd_region'])
            input_data['area'].append(val)

            assert input_data['area'][-1] != ''

        input_data['age_weights'] = [json.dumps(row.get('age_weights', [])) for row in dm['data']]  # store age_weights as json, since Pandas doesn't like arrays in arrays

        # add selected covariates
        for level in ['Country_level', 'Study_level']:
            for cv in dm['params']['covariates'][level]:
                if dm['params']['covariates'][level][cv]['rate']['value']:
                    input_data['x_%s'%cv] = [float(row.get(dismod3.utils.clean(cv), '') or 0.) for row in dm['data']]
        
        input_data = pandas.DataFrame(input_data)


        # print checks of data
        for i, row in input_data.T.iteritems():
            if pl.isnan(row['value']):
                print 'value in row %d is missing' % i
        input_data['value'] = input_data['value'].fillna(0.)

        return input_data


    @staticmethod
    def _output_template_from_gbd_json(dm):
        """ generate output template"""
        import dismod3
        output_template = {}
        for field in 'area sex year pop'.split():
            output_template[field] = []
        for level in ['Country_level', 'Study_level']:
            for cv in dm['params']['covariates'][level]:
                if dm['params']['covariates'][level][cv]['rate']['value']:
                    output_template['x_%s'%cv] = []

        for region in dismod3.settings.gbd_regions:
            for area in dm['countries_for'][dismod3.utils.clean(region)]:
                for year in dismod3.settings.gbd_years:
                    for sex in dismod3.settings.gbd_sexes:
                        sex = dismod3.utils.clean(sex)
                        output_template['area'].append(area)
                        output_template['sex'].append(sex)
                        output_template['year'].append(float(year))
                            
                        output_template['pop'].append(pl.sum(dm['population_by_age'][area, year, sex]))

                        # merge in country level covariates
                        for level in ['Country_level', 'Study_level']:
                            for cv in dm['params']['covariates'][level]:
                                if dm['params']['covariates'][level][cv]['rate']['value']:
                                    if dm['params']['covariates'][level][cv]['value']['value'] == 'Country Specific Value':
                                        if 'derived_covariate' in dm['params']:
                                            output_template['x_%s'%cv].append(dm['params']['derived_covariate'][cv].get('%s+%s+%s'%(area, year, sex)))
                                        else:
                                            output_template['x_%s'%cv].append(pl.nan)

                                    else:
                                        output_template['x_%s'%cv].append(float(dm['params']['covariates'][level][cv]['value']['value'] or 0.))
                                                
        return pandas.DataFrame(output_template)


    @staticmethod
    def _parameters_from_gbd_json(dm):
        """ copy expert priors"""
        parameters = ModelData().parameters
        old_name = dict(i='incidence', p='prevalence', rr='relative_risk', r='remission', f='excess_mortality', X='duration')
        for t in 'i p r f rr X'.split():
            parameters[t]['parameter_age_mesh'] = dm['params']['global_priors']['parameter_age_mesh']
            parameters[t]['y_maximum'] = dm['params']['global_priors']['y_maximum']
            for prior in 'smoothness heterogeneity level_value level_bounds increasing decreasing'.split():
                parameters[t][prior] = dm['params']['global_priors'][prior][old_name[t]]
        parameters['ages'] = range(dm['params']['global_priors']['parameter_age_mesh'][0], dm['params']['global_priors']['parameter_age_mesh'][-1]+1)
        return parameters


    @staticmethod
    def _hierarchy_from_gbd_json(dm):
        """ setup hierarchy and nodes_to_fit"""
        import dismod3

        superregions = [[15, 5, 9, 0, 12], [7, 8, 1], [17, 18, 19, 20], [14], [3], [4, 2, 16], [10, 11, 13, 6]]

        hierarchy = nx.DiGraph()
        nodes_to_fit = ['all']

        weight = pl.nan

        for i, superregion in enumerate(superregions):
            super_region_node = 'super-region_%d'%i
            hierarchy.add_edge('all', super_region_node, weight=weight)
            for j in superregion:
                #hierarchy.add_node(super_region_node, pop=0.)
                region_node = str(dismod3.utils.clean(dismod3.settings.gbd_regions[j]))
                nodes_to_fit.append(region_node)
                #hierarchy.add_node(region_node, pop=0.)
                hierarchy.add_edge(super_region_node, region_node, weight=weight)
                        
                for iso3 in dm['countries_for'][region_node]:
                    country_node = iso3
                    hierarchy.add_node(country_node,pop=0)
                    for year in [1990, 2005, 2010]:
                        for sex in 'male female'.split():
                            pop = sum(dm['population_by_age'][iso3, str(year), sex])
                            hierarchy.node[country_node]['pop'] += pop
                    hierarchy.add_edge(region_node, country_node, weight=weight)
                    #hierarchy.node[region_node]['pop'] += pop
                    #hierarchy.node[super_region_node]['pop'] += pop

        return hierarchy, nodes_to_fit
