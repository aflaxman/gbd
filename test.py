""" Module to test new code, and also to automatically test that working things stay working"""

# matplotlib backend setup, so that graphics can run on cluster (which doesn't have X11)
import matplotlib
matplotlib.use("AGG") 

import pylab as pl
import time

def fit_and_plot(mod, data_fname='irq_5q0.csv', image_fname='/home/j/Project/Models/space-time-smoothing/irq_test/5q0.%s.png'):
    import model
    reload(model)
    
    data = pl.csv2rec(data_fname)

    print 'generating model'
    mod_mc = eval('model.%s(data)' % mod)

    print 'fitting model with mcmc'
    mod_mc.sample(10000, 5000, 50, verbose=1)
            
    print 'summarizing results'

    import graphics
    reload(graphics)
    pl.figure(figsize=(11, 8.5), dpi=300)
    pl.clf()
    graphics.plot_prediction_over_time('IRQ', data, mod_mc.predicted, age=-1, cmap=pl.cm.GnBu, connected=False)
    graphics.plot_prediction_over_time('IRQ', data[:40], mod_mc.param_predicted, age=-1)

    #pl.plot(data.year, data.y, zorder=0,
    #        linestyle='', marker='x', mew=3, color='r', ms=8, alpha=.5)
    pl.title('IRQ')
    pl.xlabel('Time (Years)')
    pl.ylabel('$\log(_5q_0)$')
    pl.axis([1945, 2015, -1.6, .8])
    t1 = time.time()
    pl.savefig(image_fname%t1)

    try:
        print 'beta =\n%s' % mean_w_ui(mod_mc.beta)
        print 'gamma =\n%s' % mean_w_ui(mod_mc.gamma)
        print 'tau_f =\n%s' % mean_w_ui(mod_mc.tau_f)
        print 'sigma_f =\n%s' % mean_w_ui(mod_mc.sigma_f)
    except AttributeError:
        pass
    
    return mod_mc

def mean_w_ui(stoch):
    stats = stoch.stats()
    val_list = []
    for i, v in enumerate(stats['mean']):
        val_list.append('%.2f (%.2f, %.2f), ' % (v, stats['95% HPD interval'][i][0], stats['95% HPD interval'][i][1]))
    return ', '.join(val_list)
 
data_gen_models = 'fe smooth_gp_re_a'
def regenerate_data(data_model, pct=80., std=1.):
    """ Regenerate test data using specified data generation function
    Allowed models: %s """ % data_gen_models
    if data_model not in data_gen_models.split(' '):
        raise TypeError, 'Unrecognized data model "%s"; must be one of %s' % (data_model, data_gen_models)
    
    import data
    reload(data)
    eval('data.generate_%s()' % data_model)
    data.add_sampling_error(std=std)
    data.knockout_uniformly_at_random(pct=pct)

data_run_models = 'fe gp_re_a'
def evaluate_model(mod, comment='', data_fname='missing_noisy_data.csv', truth_fname='data.csv'):
    """ Run specified model on existing data (data.csv / missing_noisy_data.csv) and save results in dev_log.csv
    Existing models: %s """ % data_run_models
    if mod not in data_run_models.split(' '):
        raise TypeError, 'Unrecognized model "%s"; must be one of %s' % (mod, data_run_models)

    import model
    reload(model)

    print 'loading data'
    data = pl.csv2rec(data_fname)
    truth = pl.csv2rec(truth_fname)
    
    t0 = time.time()
    print 'generating model'
    mod_mc = eval('model.%s(data)' % mod)

    print 'fitting model with mcmc'
    mod_mc.sample(10000, 5000, 50, verbose=1)
    #mod_mc.sample(100, verbose=1)
    t1 = time.time()

    print 'summarizing results'

    import graphics
    reload(graphics)
    pl.figure(figsize=(22, 17), dpi=300)
    pl.clf()
    graphics.plot_all_predictions_over_time(data, mod_mc.predicted, more_data=truth)

    stats = mod_mc.predicted.stats()
    i_out = [i for i in range(len(data)) if pl.isnan(data.y[i])]
    rmse_abs_out = pl.rms_flat(truth.y[i_out] - stats['mean'][i_out])
    rmse_rel_out = 100*pl.rms_flat(1. - stats['mean'][i_out]/truth.y[i_out])

    i_in = [i for i in range(len(data)) if not pl.isnan(data.y[i])]
    rmse_abs_in = pl.rms_flat(truth.y[i_in] - stats['mean'][i_in])
    rmse_rel_in = 100*pl.rms_flat(1. - stats['mean'][i_in]/truth.y[i_in])

    stats = mod_mc.param_predicted.stats()
    coverage = 100*pl.sum((truth.y[i_out] >= stats['95% HPD interval'][i_out, 0]) & (truth.y[i_out] <= stats['95% HPD interval'][i_out, 1])) / float(len(i_out))

    import md5
    data_hash = md5.md5(data).hexdigest()
    data_comment = open('data.csv').readline()
    
    results = [mod, t1-t0, rmse_abs_out, rmse_rel_out, rmse_abs_in, rmse_rel_in, coverage,
               len(data), len(pl.unique(data.region)), len(pl.unique(data.country)), len(pl.unique(data.year)), len(pl.unique(data.age)), data_hash,
               t0, data_comment, comment]
    print '%s: time: %.0fs out-of-samp rmse abs=%.1f rel=%.0f in-samp rmse abs=%.1f rel=%.0f coverage=%.0f\ndata: %d rows; %d regions, %d countries %d years %d ages [data hash: %s]\n(run conducted at %f) %s - %s' % tuple(results)

    pl.savefig('/home/j/Project/Models/space-time-smoothing/images/%s.png' % t0)  # FIXME: don't hardcode path for saving images

    import csv
    f = open('dev_log.csv', 'a')
    f_csv = csv.writer(f)
    f_csv.writerow(results)
    f.close()

    return mod_mc

if __name__ == '__main__':
    import pylab as pl
    import data

    data.age_range = pl.arange(0, 81, 10)
    data.time_range = pl.arange(1980, 2005, 5)
    data.regions = pl.randint(1,22)

    time.sleep(pl.rand()*5.)
    t0 = time.time()
    data.generate_fe('test_data/%s.csv'%t0)  # included just to get good test coverage
    data.generate_smooth_gp_re_a('test_data/%s.csv'%t0, country_variation=True)

    std=5.*pl.rand(len(pl.csv2rec('test_data/%s.csv'%t0)))
    pct=25.

    print data.age_range, data.time_range, data.regions, pl.mean(std), pct
    
    data.add_sampling_error('test_data/%s.csv'%t0, 'test_data/noisy_%s.csv'%t0, std=std)
    data.knockout_uniformly_at_random('test_data/noisy_%s.csv'%t0, 'test_data/missing_noisy_%s.csv'%t0, pct=pct)

    evaluate_model('gp_re_a', 'knockout pct=%d, model that knows about varying noise and includes coutry level effects)' % pct,
                   'test_data/missing_noisy_%s.csv'%t0, 'test_data/%s.csv'%t0)
