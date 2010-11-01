""" Module to test new code, and also to automatically test that working things stay working"""

# matplotlib backend setup, so that graphics can run on cluster (which doesn't have X11)
import matplotlib
matplotlib.use("AGG") 

import pylab as pl
import time

def fit_and_plot(mod, data_fname='irq_5q0.csv', image_fname='/home/j/Project/Models/space-time-smoothing/irq_test/5q0.%s.png',
                 comment='', iter=40000):
    import model
    reload(model)
    
    data = pl.csv2rec(data_fname)

    # FIXME: this makes a big difference, but I don't understand why it would (could be prior on gp amp)
    data.x1 = (data.x1-1990.)/10.  # crude normalization of year data

    print 'generating model'
    mod_mc = eval('model.%s(data)' % mod)

    print 'fitting model with mcmc'
    mod_mc.sample(iter, iter/2, iter/2000, verbose=1)
            
    print 'summarizing results'

    import graphics
    reload(graphics)
    pl.figure()
    pl.clf()
    graphics.plot_prediction_over_time('IRQ', data, mod_mc.predicted, age=-1, cmap=pl.cm.RdYlBu, connected=False, jittered_posterior=False)
    graphics.plot_prediction_over_time('IRQ', data[:40], mod_mc.param_predicted, age=-1)

    #pl.plot(data.year, data.y, zorder=0,
    #        linestyle='', marker='x', mew=3, color='r', ms=8, alpha=.5)
    pl.title('IRQ')
    pl.xlabel('Time (Years)')
    pl.ylabel('$\log(_5q_0)$')
    pl.axis([1945, 2030, -1.8, -.5])
    pl.figtext(0, 1, '\n %s' % comment, va='top', ha='left')
    t1 = time.time()
    pl.savefig(image_fname%t1)

    try:
        for stoch in 'beta gamma sigma_f tau_f'.split(' '):
            print '%s =\n    %s\n' % (stoch, mean_w_ui(mod_mc.__getattribute__(stoch)))
    except AttributeError:
        pass
    
    return mod_mc

def mean_w_ui(stoch):
    stats = stoch.stats()
    val_list = []
    for i, v in enumerate(stats['mean']):
        val_list.append('%.2f (%.2f, %.2f) [%.2f], ' % (v, stats['95% HPD interval'][i][0], stats['95% HPD interval'][i][1], stats['mc error'][i]/abs(v)*100))
    return ', '.join(val_list)

if __name__ == '__main__':
        iter=10000
        mod_mc = fit_and_plot('gp_re_a', iter=iter,
                              comment='%dK samples, sigma_e = Exp(1)')
