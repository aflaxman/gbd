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
        self.parameters = dict(i={}, p={}, r={}, f={}, rr={}, X={}, pf={}, ages=range(101))

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

        self.input_data.to_csv(path + '/input_data.csv')
        self.output_template.to_csv(path + '/output_template.csv')
        json.dump(self.parameters, open(path + '/parameters.json', 'w'), indent=2)
        json.dump(dict(nodes=[[n, self.hierarchy.node[n]] for n in sorted(self.hierarchy.nodes())],
                       edges=[[u, v, self.hierarchy.edge[u][v]] for u,v in sorted(self.hierarchy.edges())]),
                  open(path + '/hierarchy.json', 'w'), indent=2)
        json.dump(self.nodes_to_fit, open(path + '/nodes_to_fit.json', 'w'), indent=2)

    @staticmethod
    def load(path):
        d = ModelData()

        # TODO: catch _csv.Error and retry, to give j drive time to sync
        d.input_data = pandas.DataFrame.from_csv(path + '/input_data.csv')

        # ensure that certain columns are float
        for field in 'value standard_error upper_ci lower_ci effective_sample_size'.split():
            #d.input_data.dtypes[field] = float  # TODO: figure out classy way like this, that works
            d.input_data[field] = pl.array(d.input_data[field], dtype=float)

        
        d.output_template = pandas.DataFrame.from_csv(path + '/output_template.csv')
        
        d.parameters = json.load(open(path + '/parameters.json'))

        hierarchy = json.load(open(path + '/hierarchy.json'))
        d.hierarchy.add_nodes_from(hierarchy['nodes'])
        d.hierarchy.add_edges_from(hierarchy['edges'])

        d.nodes_to_fit = json.load(open(path + '/nodes_to_fit.json'))

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

        # remove any rows with 'ignore' columns set to 1
        dm['data'] = [d for d in dm['data'] if not (d.get('Ignore') or d.get('ignore'))]

        # remove any data with type-specific heterogeneity set to Unusable
        if 'global_priors' in dm['params']:
            for t in dm['params']['global_priors']['heterogeneity']:
                if dm['params']['global_priors']['heterogeneity'][t] == 'Unusable':
                    print '%s has heterogeneity unusable, dropping %d rows' % (t, len([d for d in dm['data'] if d['data_type'] == t + ' data']))
                    dm['data'] = [d for d in dm['data'] if d['data_type'] != t + ' data']

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
                         'prevalence x excess-mortality data': 'pf', 'all-cause mortality data': 'm_all', 'relative-risk data': 'rr',
                         'duration data': 'X', 'smr data': 'smr', 'cause-specific mortality data': 'csmr', 'mortality data': 'm_with'}
        input_data['data_type'] = [new_type_name[row['data_type']] for row in dm['data']]

        for field in 'value standard_error lower_ci upper_ci'.split():
            input_data[field] = []
            for row in dm['data']:
                val = row.get(field, '')
                if val == '':
                    val = pl.nan
                else:
                    val = float(val) / float(row.get('units', '1').replace(',', ''))
                input_data[field].append(val)

        input_data['area'] = []
        for row in dm['data']:
            val = row.get('country_iso3_code', '')
            if val == '' or val == 'all':
                val = dismod3.utils.clean(row['gbd_region'])
            input_data['area'].append(val)

            assert input_data['area'][-1] != ''

        input_data['age_weights'] = [';'.join(['%.4f'%w for w in row.get('age_weights', [])]) for row in dm['data']]  # store age_weights as semi-colon delimited text, since Pandas doesn't like arrays in arrays and doesn't save comma-separated fields correctly

        # add selected covariates
        if 'covariates' in dm['params']:
            for level in ['Country_level', 'Study_level']:
                for cv in dm['params']['covariates'][level]:
                    if dm['params']['covariates'][level][cv]['rate']['value']:
                        input_data['x_%s'%cv] = [float(row.get(dismod3.utils.clean(cv), '') or 0.) for row in dm['data']]

                    # also include column of input data for 'z_%s'%cv if it is requested
                    if dm['params']['covariates'][level][cv]['error']['value']:
                        input_data['z_%s'%cv] = [float(row.get(dismod3.utils.clean(cv), '') or 0.) for row in dm['data']]

        input_data = pandas.DataFrame(input_data)


        # replace age_end 1 with age_end 0, correcting a common mistake in data entry
        i = (input_data['age_start']==0) & (input_data['age_end']==1)
        if i.sum() > 0:
            print 'WARNING: correcting age_end in %d rows that have age_start == 0, age_end == 1 (old format uses "demographic" notation)' % i.sum()
            input_data['age_end'][i] = 0

        # replace triple underscores with single underscore, a problem with consistency in the spacing in "North Africa / Middle East"
        input_data['area'] = [a.replace('___', '_') for a in input_data['area']]

        # print checks of data
        for i, row in input_data.T.iteritems():
            if pl.isnan(row['value']):
                print 'WARNING: value in row %d is missing' % i
        input_data = input_data[~pl.isnan(input_data['value'])]

        return input_data


    @staticmethod
    def _output_template_from_gbd_json(dm):
        """ generate output template"""
        import dismod3
        output_template = {}
        for field in 'area sex year pop'.split():
            output_template[field] = []
        if 'covariates' in dm['params']:
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
                        if 'covariates' in dm['params']:
                            for level in ['Country_level', 'Study_level']:
                                for cv in dm['params']['covariates'][level]:
                                    if dm['params']['covariates'][level][cv]['rate']['value']:
                                        
                                        if level == 'Country_level' and dm['params']['covariates'][level][cv]['value']['value'] == '':
                                            # people usually mean CSV, so interpret blanks to mean this
                                            dm['params']['covariates'][level][cv]['value']['value'] = 'Country Specific Value'
                                            
                                        if dm['params']['covariates'][level][cv]['value']['value'] == 'Country Specific Value':
                                            if 'derived_covariate' in dm['params'] and cv in dm['params']['derived_covariate']:
                                                output_template['x_%s'%cv].append(dm['params']['derived_covariate'][cv].get('%s+%s+%s'%(area, year, sex)))
                                                
                                            else:
                                                raise KeyError, 'covariate %s not found for output template (did you set a reference value? did you "Calculate covariates for model data"?)' % cv

                                        else:
                                            output_template['x_%s'%cv].append(float(dm['params']['covariates'][level][cv]['value']['value'] or 0.))
                                                
        return pandas.DataFrame(output_template)


    @staticmethod
    def _parameters_from_gbd_json(dm):
        """ copy expert priors"""
        parameters = ModelData().parameters
        old_name = dict(i='incidence', p='prevalence', rr='relative_risk', r='remission', f='excess_mortality', X='duration', pf='prevalence_x_excess-mortality')
        for t in 'i p r f rr X pf'.split():
            if 'global_priors' in dm['params']:
                parameters[t]['parameter_age_mesh'] = dm['params']['global_priors']['parameter_age_mesh']
                parameters[t]['y_maximum'] = dm['params']['global_priors']['y_maximum']
                for prior in 'smoothness heterogeneity level_value level_bounds increasing decreasing'.split():
                    if old_name[t] in dm['params']['global_priors'][prior]:
                        parameters[t][prior] = dm['params']['global_priors'][prior][old_name[t]]
                    elif old_name[t] == 'prevalence_x_excess-mortality':
                        parameters[t][prior] = dm['params']['global_priors'][prior]['excess_mortality']

                # make 1000 effectively infinite, because the gui only goes up to 1000
                if parameters[t]['level_bounds']['upper'] == 1000.:
                    parameters[t]['level_bounds']['upper'] = 1e6
            parameters[t]['fixed_effects'] = {}
            parameters[t]['random_effects'] = {}

        if 'global_priors' in dm['params']:
            parameters['ages'] = range(dm['params']['global_priors']['parameter_age_mesh'][0], dm['params']['global_priors']['parameter_age_mesh'][-1]+1)

        for t in 'i p r f'.split():
            key = 'sex_effect_%s' % old_name[t]
            if key in dm['params']:
                prior = dm['params'][key]
                parameters[t]['fixed_effects']['x_sex'] = dict(dist='normal', mu=pl.log(prior['mean']),
                                                               sigma=(pl.log(prior['upper_ci']) - pl.log(prior['lower_ci']))/(2*1.96))
            key = 'region_effect_%s' % old_name[t]
            if key in dm['params']:
                prior = dm['params'][key]
                for iso3 in dm['countries_for']['World']:
                    parameters[t]['random_effects'][iso3] = dict(dist='TruncatedNormal', mu=0., sigma=prior['std'], lower=-2*prior['std'], upper=2*prior['std'])

            # include alternative prior on sigma_alpha based on heterogeneity
            for i in range(5):  # max depth of hierarchy is 5
                effect = 'sigma_alpha_%s_%d'%(t,i)
                #parameters[t]['random_effects'][effect] = dict(dist='TruncatedNormal', mu=.01, sigma=.01, lower=.01, upper=.05)
                #if 'heterogeneity' in parameters[t]:
                #    if parameters[t]['heterogeneity'] == 'Moderately':
                #        parameters[t]['random_effects'][effect] = dict(dist='TruncatedNormal', mu=.05, sigma=.05, lower=.01, upper=1.)
                #    elif parameters[t]['heterogeneity'] == 'Very':
                #        parameters[t]['random_effects'][effect] = dict(dist='TruncatedNormal', mu=.01, sigma=.01, lower=.002, upper=.2)
            
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
