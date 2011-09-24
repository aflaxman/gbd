""" Data Handling Class for DisMod III"""

import pandas
import networkx as nx

class Data:
    def __init__(self):
        self.input_data = pandas.DataFrame(columns='data_type value area sex age_start age_end year_start year_end standard_error effective_sample_size lower_ci upper_ci age_weights'.split())
        self.output_template = pandas.DataFrame()
        self.parameters = dict()
        self.areas_hierarchy = nx.DiGraph()
        self.areas_to_fit = list()
        
