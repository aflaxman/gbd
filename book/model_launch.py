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

draws = sys.argv[1]
area = 'europe_western'
data_type = 'p'

m = 29561
os.system('/usr/local/bin/SGE/bin/lx24-amd64/qsub -cwd /homes/peterhm/gbd/book/model_comparison.sh %s %d' %(draws, m))
# # load best models spread sheet
# bm_path = '../../GBD/dalynator/yld/best_models.csv'
# bm_csv = pandas.read_csv(bm_path,index_col=None)

# dismod_models = bm_csv.groupby('dismod_model_number').apply(lambda df: df.ix[df.index[0], 'outcome_name'])
# dismod_models = dismod_models.drop([0], axis=0)

# model_list = []
# for i, m in enumerate(dismod_models.index):
    # m = int(m)
    # try:
        # model = mu.load_new_model(m, area, data_type)
        # if len(model.input_data['data_type'].index) >= 100: model_list.append(m)
        # # submit shell
        # os.system('/usr/local/bin/SGE/bin/lx24-amd64/qsub /homes/peterhm/gbd/book/model_comparison.sh %s %d' %(draws, m))
    # except IOError:
        # pass



