import matplotlib
matplotlib.use("AGG")

import optparse
import dismod3
from pylab import *

def plot_x_eq_y():
    l,r,b,t = axis()
    m = max(r,t)
    n = min(l,b)
    plot([n, m], [n, m], '--', color='grey')
    axis([l,r,b,t])

def plot_comparison(dm_old, dm_new, type):
    for k in dm_new.params.get('mcmc_mean', {}).keys():
        if k.startswith(type):
            x=dm_old.get_mcmc('mean', k)
            y=dm_new.get_mcmc('mean', k)
            if len(x) != len(y):
                continue
            if rms_flat(array(x) - array(y)) < 1.e-7:
                continue
            plot(x, y, '.', label=k)
    title('%d/%d - %s' % (dm_old.id, dm_new.id, dm_old.get_condition()))
    xlabel('old %s'%type)
    ylabel('new %s'%type)
    #legend(loc=(1,0))
    plot_x_eq_y()
    grid()


import mare

rmse_dict = {}

def compare_residuals(old_id, new_id):
    df_new = mare.upload_fits.merge_data_csvs(new_id)
    df_old = mare.upload_fits.merge_data_csvs(old_id)
    
    df = df_new.join(df_old,lsuffix='_new', rsuffix='_old')

    plot(array(df['residual_old']), array(df['residual_new']), 'ks', alpha=.25)

    title('%d/%d - %s' % (old_id, new_id, 'residuals'))
    xlabel('old residual')
    ylabel('new residual')
    
    rmse_str = 'MAE:\n'
    rmse_str += '  old: %.6f\n'%median(df['abs_residual_old'])
    rmse_str += '  new: %.6f\n'%median(df['abs_residual_new'])
    
    rmse_dict[id, 'old'] = median(df['abs_residual_old'])
    rmse_dict[id, 'new'] = median(df['abs_residual_new'])
    
    l,r,b,t = axis()
    text(r,b,rmse_str, ha='right', va='bottom')
    #legend(loc=(1,0))
    plot_x_eq_y()
    grid()
    
    return df

def plot_age_pattern(dm_old, dm_new, type, ax_old, ax_new):
    for k in dm_new.params.get('mcmc_mean', {}).keys():
        if k.startswith(type):
            y_old=dm_old.get_mcmc('mean', k)
            y_new=dm_new.get_mcmc('mean', k)
            if len(y_new) != 101 or len(y_old) != 101:
                continue
            if rms_flat(array(y_old) - array(y_new)) < 1.e-7:
                continue
            ax_old.plot(range(101), y_old, '-', label=k)
            ax_new.plot(range(101), y_new, '-', label=k)
    ax_old.set_title('old')
    ax_new.set_title('new')
    for ax in [ax_old, ax_new]:
        ax.set_xlabel('age')
        ax.set_ylabel(type)
        ax.grid()
        
if __name__ == '__main__':
    usage = 'usage: %prog [options] old_id new_id'
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.error('incorrect number of arguments')

    try:
        old_id = int(args[0])
        new_id = int(args[1])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    dm_old = dismod3.disease_json.load_disease_model(old_id)
    dm_new = dismod3.disease_json.load_disease_model(new_id)

    type = 'prevalence'
    figure(figsize=(22,12), dpi=600)
    clf()
    subplot(2,2,1)
    plot_comparison(dm_old, dm_new, type)
    
    subplot(2,2,3)
    try:
        df = compare_residuals(old_id, new_id)
        
    except Exception, e:
        print e
        text(.5, .5, 'residuals failed or not yet ready')
        
        
    plot_age_pattern(dm_old, dm_new, type, subplot(2,2,2), subplot(2,2,4))

    subplots_adjust(hspace=.5)
    savefig('/home/j/Project/dismod/comparison_%d_%d.png'%(old_id, new_id))
