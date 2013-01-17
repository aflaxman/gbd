import sys
sys.path += ['../gbd', '../gbd/book', '../dm3-computation_only/', '../dm3-computation_only/book']
import pylab as pl
import pymc as mc
import pandas

import dismod3
reload(dismod3)

import book_graphics
reload(book_graphics)
import matplotlib as mpl

# set font
book_graphics.set_font()

def load_new_model():
    orig_model = dismod3.data.load('/home/j/Project/dismod/notebooks/hep_c')
    return orig_model

def my_axis(ymax):
    pl.axis([-5,105,-ymax/10.,ymax])

best_model = load_new_model()

# finding hierarchy moving downstream (super-region -> region -> country)
sr2r = {}
r2c = {}

sr = best_model.hierarchy['all'].keys()
for sr_i in range(len(sr)):
    sr2r[sr[sr_i]] = best_model.hierarchy[sr[sr_i]].keys() # sr <- r
    
    r = best_model.hierarchy[sr[sr_i]].keys()
    for r_i in range(len(r)):
        r2c[r[r_i]] = best_model.hierarchy[r[r_i]].keys() # r <- c
    
# finding hierarchy moving upstream (country -> region -> super-region)
c2r = {}
r2sr = {}

sr = best_model.hierarchy['all'].keys()
for sr_i in range(len(sr)):
    r = best_model.hierarchy[sr[sr_i]].keys()
    for r_i in range(len(r)):
        r2sr[r[r_i]] = sr[sr_i] # r <- sr 
        c = best_model.hierarchy[r[r_i]].keys()
        for c_i in range(len(c)):
            c2r[c[c_i]] = r[r_i] # c <- r

def levels(area):
    levels = {}
    if len(best_model.hierarchy[area].keys()) != 0: # if there are no keys, the area is not a country
        try: # if area is a region
            levels['countries'] = r2c[area]
            levels['region'] = area
            levels['superregion'] = r2sr[area]
        except KeyError: # if area a super-region
            levels['superregion'] = area
            levels['region'] = sr2r[area]
    else: # if area is a country
        levels['country'] = area
        levels['region'] = c2r[area]
        levels['superregion'] = r2sr[levels['region']]
        
    return levels            
    
# figure hepc-EGY_v_JOR
egy_model = load_new_model()
egy_model.keep(areas=['EGY'])

jor_model = load_new_model()
jor_model.keep(areas=['JOR'])

pl.figure(**book_graphics.half_page_params)

pl.subplot(1,2,1)
dismod3.graphics.plot_data_bars(egy_model.get_data('p'))
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .15, .30, .45, .60],[0, 15, 30, 45, 60])
my_axis(.7)
book_graphics.subtitle('(a)')


pl.subplot(1,2,2)
dismod3.graphics.plot_data_bars(jor_model.get_data('p'))
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .0025, .005, .0075, .01], [0, .25, .5, .75, 1])
my_axis(.012)
book_graphics.subtitle('(b)')


pl.subplots_adjust(wspace=.35, top=.99, bottom=.14)

pl.savefig('book/graphics/hepc-EGY_v_JOR.pdf')
pl.savefig('book/graphics/hepc-EGY_v_JOR.png')

# figure hepc-region_v_EGY_v_JOR
regional = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_hepc_region.csv')
pl.figure(**book_graphics.full_page_params)
     
pl.plot(pl.array(regional['NAME']), 'k-', linewidth=2, label='Region')
pl.plot(pl.array(regional['EGY']), 'k--', linewidth=2, label='Egypt')
pl.plot(pl.array(regional['JOR']), 'k:', linewidth=2, label='Jordon')

pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .07, .14, .21, .28], [0, 7, 14, 21, 28])
my_axis(.3)
pl.legend(loc='upper right', fancybox=True, shadow=True)
pl.subplots_adjust(hspace=.35)
pl.subplots_adjust(wspace=.35)


pl.savefig('book/graphics/hepc-region_v_EGY_v_JOR.pdf')
pl.savefig('book/graphics/hepc-region_v_EGY_v_JOR.png')

# figure hepc-tree_plot_global_hetero
data = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_hepc.csv')

model1 = load_new_model()
model1.vars += dismod3.ism.age_specific_rate(model1, 'p')

captions = dict(EGY = 'Egypt', 
                JOR = 'Jordan',
                SAU = 'Saudi Arabia',
                IRQ = 'Iraq',
                IRN = 'Iran',
                YEM = 'Yemen',
                TUR = 'Turkey',
                SYR = 'Syria',
                TUN = 'Tunisia',
                north_africa_middle_east = 'North Africa and Middle East',
                superregion_6 = 'Latin America and Caribbean',
                superregion_5 = 'Southeast Asia, East Asia, and Oceania',
                superregion_4 = 'South Asia',
                superregion_3 = 'North Africa and Middle East',
                superregion_2 = 'Sub-Saharan Afica',
                superregion_1 = 'Central Europe, Eastern Europe, and Central Asia',
                superregion_0 = 'High-income'
                )
                
pl.figure(**book_graphics.full_page_params)

import networkx as nx
rate_type='p'
area='EGY'

hierarchy = model1.hierarchy
covariate = 'U'
effect = 'alpha'

#pl.figure(figsize=(11, 12))
l=-2.5
r=3

# list of all covariate names
cov_name = list(model1.vars['p'][covariate].columns)
# list of countries, regions, and super-regions of interest
hlevels = levels(area)
hlvl_list = []
hlvl_list.extend(sr2r.keys())
hlvl_list.extend(sr2r[hlevels['superregion']])
if hlevels['superregion'] != area: hlvl_list.extend(r2c[hlevels['region']])
# intersection of covariates and areas of interest
plot_list = list(set(cov_name) & set(hlvl_list))

if isinstance(model1.vars['p'].get(effect), list):
    index = sorted(pl.arange(len(cov_name)),
                   key=lambda i: str(cov_name[i] in hierarchy and nx.shortest_path(hierarchy, 'all', cov_name[i]) or cov_name[i]))

# index of area of interest in covariate list
plot_list_ix = []
for i in range(len(plot_list)): 
    plot_list_ix.append(cov_name.index(plot_list[i]))

# correctly ordering geographic regions
plot_list_ix_order = []
for i in range(len(index)):
    if (index[i] in plot_list_ix) == 1: plot_list_ix_order.append(index[i]) # if true, append
    
for y, i in enumerate(plot_list_ix_order):
    x1 = data.ix[cov_name[i],'moderately']
    x2 = data.ix[cov_name[i],'slightly'] 
    x3 = data.ix[cov_name[i],'very']
    
    xerr1 = pl.array([x1 - pl.atleast_2d(data.ix[cov_name[i],'moderately_l']),
                      pl.atleast_2d(data.ix[cov_name[i],'moderately_u']) - x1])
    xerr2 = pl.array([x2 - pl.atleast_2d(data.ix[cov_name[i],'slightly_l']),
                      pl.atleast_2d(data.ix[cov_name[i],'slightly_u']) - x2])
    xerr3 = pl.array([x3 - pl.atleast_2d(data.ix[cov_name[i],'very_l']),
                      pl.atleast_2d(data.ix[cov_name[i],'very_u']) - x3])
        
    pl.errorbar(x2, 2.75*y+.45, xerr=xerr2, fmt='ko', mec='w')#, label = 'Slightly'
    pl.errorbar(x1, 2.75*y, xerr=xerr1, fmt='k^', mec='w')#, label = 'Moderately'
    pl.errorbar(x3, 2.75*y-.45, xerr=xerr3, fmt='ks', mec='w')#, label = 'Very'
        
pl.plot([-10],[-10], 'ko-', mec='w', label = '$\delta \sim \\mathrm{Uniform}(9,$ $81)$')
pl.plot([-10],[-10], 'k^-', mec='w', label = '$\delta \sim \\mathrm{Uniform}(3,$ $27)$')
pl.plot([-10],[-10], 'ks-', mec='w', label = '$\delta \sim \mathrm{Uniform}(1,$ $9)$')
l,r,b,t = pl.axis()
r = 3.0
l = -3.5
pl.vlines([0], b-.5, t+.5)
pl.xticks([l, 0, r])
pl.yticks([])
   
# adds names to data    
for y, i in enumerate(plot_list_ix_order):
    spaces = cov_name[i] in hierarchy and len(nx.shortest_path(hierarchy, 'all', cov_name[i])) or 0
    # removes - from cov_name lists and gives proper caption name
    w = list(cov_name[i])
    for c in range(len(w)):
        if w[c] == '-': w[c] = ''
    cov_name[i] = captions["".join(w)]
    pl.text(l, 2.75*y, '%s%s' % ('  '*spaces, cov_name[i]), va='center', ha='left')

pl.axis([l, r, -2.75, t-2.75])
pl.xticks([-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5])

pl.legend(title='Prior on the overdispersion, $\delta$',numpoints=1, fancybox=True, shadow=True)

pl.savefig('book/graphics/hepc-tree_plot_global_hetero.pdf')
pl.savefig('book/graphics/hepc-tree_plot_global_hetero.png')


pl.show()