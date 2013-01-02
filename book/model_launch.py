import os
import sys
sys.path += ['.', '..', '/homes/peterhm/gbd/', '/homes/peterhm/gbd/book'] 
import pylab as pl
import pymc as mc
import pandas
import random

import dismod3
reload(dismod3)

import model_utilities as mu
reload(mu)

replicates = int(sys.argv[1])
area = 'europe_western'
data_type = 'p'

rate_types = ['neg_binom'] #['log_offset', 'normal', 'log_normal', 'binom']

# load best models spread sheet
bm_path = '/snfs1/Project/GBD/dalynator/yld/best_models.csv'
bm_csv = pandas.read_csv(bm_path,index_col=None)

dismod_models = bm_csv.groupby('dismod_model_number').apply(lambda df: df.ix[df.index[0], 'outcome_name'])
dismod_models = dismod_models.drop([0], axis=0)

model_list = []
name_list = []
for m in [38393, 38395, 38398]: #dismod_models.index:
    m = int(m)
    try:
        # check that model has more than 100 prevalence points 
        model = mu.load_new_model(m, area, data_type)
        if len(model.input_data['data_type'].index) >= 100: model_list.append(m)
        for r in rate_types:
            for i in range(replicates):
                # create unique id for hold
                name = r + str(m) + str(i)
                name_list.append(name)
                # submit shell
                os.system('/usr/local/bin/SGE/bin/lx24-amd64/qsub -cwd -N ' + name + ' /homes/peterhm/gbd/book/model_comparison.sh %d %s %d' %(m, r, i))
    except IOError:
        pass

# save information about rate types and models
model_info = pandas.DataFrame(model_list, columns=['model_list'])
model_info.to_csv('/homes/peterhm/gbd/book/validity/model_list.csv')
model_types = pandas.DataFrame(rate_types, columns=['rate_types'])
model_types.to_csv('/homes/peterhm/gbd/book/validity/model_types.csv')
    
# joining all jobs into files
hold_str = '-hold_jid %s ' % ','.join(name_list)
os.system('/usr/local/bin/SGE/bin/lx24-amd64/qsub -cwd ' + hold_str + '/homes/peterhm/gbd/book/model_join.sh %d' %(replicates))




