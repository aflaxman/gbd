# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pandas 

# <codecell>

# load data
try:
    data = pandas.read_csv('/snfs1/DATA/IHME_COUNTRY_CODES/IHME_COUNTRYCODES.CSV', index_col=None)
except:
    data = pandas.read_csv('J:/DATA/IHME_COUNTRY_CODES/IHME_COUNTRYCODES.CSV', index_col=None)

# <codecell>

# keep superregion, region, country and iso3 code
data = data.ix[:,['ihme_indic_country', 'gbd_super_region_name', 'gbd_region', 'countryname_ihme', 'iso3']]

# <codecell>

# keep only countries for which IHME reports
data = data[data.ix[:,'ihme_indic_country'] == 1]

# drop duplicates
data = data.drop_duplicates('iso3')

# drop IHME indicator
data = data.drop('ihme_indic_country', axis=1)

# <codecell>

# IHME writing style guidelines
# superregions
data.ix[data['gbd_super_region_name'] == 'East Asia/Pacific', 'gbd_super_region_name'] = 'Southeast Asia, East Asia, and Oceania'
data.ix[data['gbd_super_region_name'] == 'Eastern Europe/Central Asia', 'gbd_super_region_name'] = 'Central Europe, Eastern Europe, and Central Asia'
data.ix[data['gbd_super_region_name'] == 'High Income', 'gbd_super_region_name'] = 'High-income'
data.ix[data['gbd_super_region_name'] == 'Latin America/Caribbean', 'gbd_super_region_name'] = 'Latin America and Caribbean'
data.ix[data['gbd_super_region_name'] == 'North Africa/Middle East', 'gbd_super_region_name'] = 'North Africa and Middle East'

# regions
data.ix[data['gbd_region'] == 'Asia Pacific, High Income', 'gbd_region'] = 'High-income Asia Pacific'
data.ix[data['gbd_region'] == 'Europe, Western', 'gbd_region'] = 'Western Europe'
data.ix[data['gbd_region'] == 'North America, High Income', 'gbd_region'] = 'High-income North America'
data.ix[data['gbd_region'] == 'Europe, Central', 'gbd_region'] = 'Central Europe'
data.ix[data['gbd_region'] == 'Latin America, Southern', 'gbd_region'] = 'Southern Latin America'
data.ix[data['gbd_region'] == 'Europe, Eastern', 'gbd_region'] = 'Eastern Europe'
data.ix[data['gbd_region'] == 'Asia, East', 'gbd_region'] = 'East Asia'
data.ix[data['gbd_region'] == 'Latin America, Tropical', 'gbd_region'] = 'Tropical Latin America'
data.ix[data['gbd_region'] == 'Latin America, Central', 'gbd_region'] = 'Central Latin America'
data.ix[data['gbd_region'] == 'Asia, Southeast', 'gbd_region'] = 'Southeast Asia'
data.ix[data['gbd_region'] == 'Asia, Central', 'gbd_region'] = 'Central Asia'
data.ix[data['gbd_region'] == 'Latin America, Andean', 'gbd_region'] = 'Andean Latin America'
data.ix[data['gbd_region'] == 'North Africa/Middle East', 'gbd_super_region_name'] = 'North Africa and Middle East'
data.ix[data['gbd_region'] == 'Asia, South', 'gbd_region'] = 'South Asia'
data.ix[data['gbd_region'] == 'Sub-Saharan Africa, Southern', 'gbd_region'] = 'Southern sub-Saharan Africa'
data.ix[data['gbd_region'] == 'Sub-Saharan Africa, East', 'gbd_region'] = 'Eastern sub-Saharan Africa'
data.ix[data['gbd_region'] == 'Sub-Saharan Africa, Central', 'gbd_region'] = 'Central sub-Saharan Africa'
data.ix[data['gbd_region'] == 'Sub-Saharan Africa, West', 'gbd_region'] = 'Western sub-Saharan Africa'

# countries
data.ix[data['iso3'] == 'BRN', 'countryname_ihme'] = 'Brunei'
data.ix[data['iso3'] == 'CIV', 'countryname_ihme'] = 'C\\^ote d\'Ivoire'
data.ix[data['iso3'] == 'COD', 'countryname_ihme'] = 'Democratic Republic of Congo'
data.ix[data['iso3'] == 'FSM', 'countryname_ihme'] = 'Federated States of Micronesia'
data.ix[data['iso3'] == 'GBR', 'countryname_ihme'] = 'UK'
data.ix[data['countryname_ihme'] == 'Iran, Islamic Republic of', 'countryname_ihme'] = 'Iran'
data.ix[data['countryname_ihme'] == 'Korea, Republic of', 'countryname_ihme'] = 'South Korea'
data.ix[data['iso3'] == 'LAO', 'countryname_ihme'] = 'Laos'
data.ix[data['countryname_ihme'] == 'Libyan Arab Jamahiriya', 'countryname_ihme'] = 'Libya'
data.ix[data['countryname_ihme'] == 'Macedonia, the Former Yugoslav Republic of', 'countryname_ihme'] = 'Macedonia'
data.ix[data['iso3'] == 'MMR', 'countryname_ihme'] = 'Burma'
data.ix[data['iso3'] == 'PRK', 'countryname_ihme'] = 'North Korea'
data.ix[data['iso3'] == 'PSE', 'countryname_ihme'] = 'Occupied Palestinian Territory'
data.ix[data['countryname_ihme'] == 'Russian Federation', 'countryname_ihme'] = 'Russia'
data.ix[data['iso3'] == 'STP', 'countryname_ihme'] = 'S\\~ao Tom\\\'e and Pr\\\'ncipe'
data.ix[data['countryname_ihme'] == 'Syrian Arab Republic', 'countryname_ihme'] = 'Syria'
data.ix[data['countryname_ihme'] == 'Taiwan, Province of China', 'countryname_ihme'] = 'Taiwan'
data.ix[data['countryname_ihme'] == 'Tanzania, United Republic of', 'countryname_ihme'] = 'Tanzania'
data.ix[data['iso3'] == 'USA', 'countryname_ihme'] = 'USA'
data.ix[data['countryname_ihme'] == 'Viet Nam', 'countryname_ihme'] = 'Vietnam'


# <codecell>

# create latex table
print data.sort(['gbd_super_region_name', 'gbd_region']).to_latex(index=False)

# <codecell>


