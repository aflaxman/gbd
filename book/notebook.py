# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pandas 

# <codecell>

# load data
data = pandas.read_csv('/snfs1/DATA/IHME_COUNTRY_CODES/IHME_COUNTRYCODES.CSV')

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
data.ix[data['gbd_super_region_name'] == 'East Asia/Pacific', 'gbd_super_region_name'] = 'Southeast Asia, East Asia, and Oceania'
data.ix[data['gbd_super_region_name'] == 'Eastern Europe/Central Asia', 'gbd_super_region_name'] = 'Central Europe, Eastern Europe, and Central Asia'
data.ix[data['gbd_super_region_name'] == 'Latin America/Caribbean', 'gbd_super_region_name'] = 'Latin America and Caribbean'
data.ix[data['gbd_super_region_name'] == 'North Africa/Middle East', 'gbd_super_region_name'] = 'North Africa and Middle East'
data.ix[data['gbd_super_region_name'] == 'High Income', 'gbd_super_region_name'] = 'High-income'

data.ix[data['gbd_region'] == 'Europe, Eastern', 'gbd_region'] = 'Eastern Europe'
data.ix[data['gbd_region'] == 'Europe, Western', 'gbd_region'] = 'Western Europe'
data.ix[data['gbd_region'] == 'Latin America, Southern', 'gbd_region'] = 'Latin America, South'
data.ix[data['gbd_region'] == 'Sub-Saharan Africa, Southern', 'gbd_region'] = 'Sub-Saharan Africa, South'

data.ix[data['countryname_ihme'] == 'Iran, Islamic Republic of', 'countryname_ihme'] = 'Iran'
data.ix[data['countryname_ihme'] == 'Libyan Arab Jamahiriya', 'countryname_ihme'] = 'Libya'
data.ix[data['countryname_ihme'] == 'Macedonia, the Former Yugoslav Republic of', 'countryname_ihme'] = 'Macedonia'
data.ix[data['countryname_ihme'] == 'Russian Federation', 'countryname_ihme'] = 'Russia'
data.ix[data['countryname_ihme'] == 'Korea, Republic of', 'countryname_ihme'] = 'South Korea'
data.ix[data['countryname_ihme'] == 'Syrian Arab Republic', 'countryname_ihme'] = 'Syria'
data.ix[data['countryname_ihme'] == 'Taiwan, Province of China', 'countryname_ihme'] = 'Taiwan'
data.ix[data['countryname_ihme'] == 'Tanzania, United Republic of', 'countryname_ihme'] = 'Tanzania'
data.ix[data['countryname_ihme'] == 'Viet Nam', 'countryname_ihme'] = 'Vietnam'
data.ix[data['iso3'] == 'PRK', 'countryname_ihme'] = 'North Korea'
data.ix[data['iso3'] == 'LAO', 'countryname_ihme'] = 'Laos'

# <codecell>

# create latex table
print data.sort(['gbd_super_region_name', 'gbd_region']).to_latex(index=False)

# <codecell>


