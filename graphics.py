import pylab as pl
import pymc as mc
import pandas
import networkx as nx


def plot_data_bars(df):
    """ Plot some data bars
    Input
    -----
    df : pandas.DataFrame with columns age_start, age_end, value
    """
    data_bars = zip(df['age_start'], df['age_end'], df['value'])

    # show at most 500 bars, to keep things fast
    if len(data_bars) > 500:
        import random
        data_bars = random.sample(data_bars, 500)

    # make lists of x and y points, faster than ploting each bar
    # individually
    x = []
    y = []
    for a_0i, a_1i, p_i in data_bars:
        x += [a_0i, a_1i, pl.nan]
        y += [p_i, p_i, pl.nan]

    pl.plot(x, y, 'ks-', mew=1, mec='w', ms=4)
    
def plot_fit(model, vars, emp_priors, posteriors):
    """ plot results of a fit"""
    pl.figure()
    ages = vars['i']['ages']  # not all data models have an ages key, but incidence always does
    for j, t in enumerate('i r f p rr pf'.split()):
        pl.subplot(2, 3, j+1)
        plot_data_bars(model.input_data[model.input_data['data_type'] == t])
        try:
            pl.plot(ages, vars[t]['mu_age'].stats()['mean'], 'k-', linewidth=2)
            pl.plot(ages, vars[t]['mu_age'].stats()['95% HPD interval'], 'k--')
        except TypeError:
            print 'Could not generate output statistics'
        pl.plot(ages, posteriors[t], color='b', linewidth=1)
        if t in emp_priors:
            pl.plot(ages, emp_priors[t], color='r', linewidth=1)
        pl.title(t)

def plot_one_type(model, vars, emp_priors, t):
    """ plot results of fit for one data type only"""
    pl.figure()
    plot_data_bars(model.input_data[model.input_data['data_type'] == t])
    pl.plot(vars['ages'], vars['mu_age'].stats()['mean'], 'k-', linewidth=2)
    pl.plot(vars['ages'], vars['mu_age'].stats()['95% HPD interval'], 'k--')
    pl.plot(vars['ages'], emp_priors[t], color='r', linewidth=1)
    pl.title(t)

def plot_one_ppc(vars, t):
    """ plot data and posterior predictive check"""
    pl.figure()
    pl.title(t)

    p = vars['p_obs'].value.__array__()
    n = vars['p_obs'].parents['n'].__array__()
    i = (-p).argsort()

    y = pl.arange(len(i), dtype=float)
    x = p[i]
    xerr = pl.sqrt(p[i] * (1-p[i]) / n[i])
    #pl.errorbar(x, y, xerr=xerr, fmt='ks', mec='w', label='Observed Data')
    pl.plot(x, y, 'ks', mec='w', label='Observed Data')

    y += .2
    stats = vars['p_pred'].stats()
    x = stats['mean'][i]
    xerr = [x - pl.atleast_2d(stats['95% HPD interval'])[i,0],
            pl.atleast_2d(stats['95% HPD interval'])[i,1] - x]
    pl.errorbar(x, y, xerr=xerr, fmt='ko', mec='w', label='Predicted Data')

    pl.yticks([])
    pl.legend(fancybox=True, shadow=True)

    l,r,b,t = pl.axis()
    pl.axis([0, r, b-.5, t+.5])



def plot_effects(vars):
    """ plot the effect coefficients for a consistent fit"""
    pl.figure()

    # count how many data models have effect coefficients
    rows = 0
    for type in 'i r f p rr pf'.split():
        if isinstance(vars[type].get('beta'), mc.Stochastic):
            rows += 1

    tile = 1
    for type in 'i r f p rr pf'.split():
        for i, (covariate, effect) in enumerate([['U', 'alpha'], ['X', 'beta']]):
            if isinstance(vars[type].get(effect), mc.Stochastic):
                pl.subplot(rows, 2, tile)
                pl.title('%s_%s' % (effect, type))

                stats = vars[type][effect].stats()
                if stats:
                    x = pl.atleast_1d(stats['mean'])
                    y = range(len(x))

                    xerr = [x - pl.atleast_2d(stats['95% HPD interval'])[:,0],
                            pl.atleast_2d(stats['95% HPD interval'])[:,1] - x]
                    pl.errorbar(x, y, xerr=xerr, fmt='bs', mec='w')

                    l,r,b,t = pl.axis()
                    pl.vlines([0], b-.5, t+.5)
                    pl.hlines(y, l, r, linestyle='dotted')
                    pl.axis([l,r,b-.5, t+.5])
                    pl.xticks([l, 0, r])
                    pl.yticks([])
                    for y_i, cov_i in enumerate(list(vars[type][covariate].columns)):
                        pl.text(l, y_i, ' %s\n' % cov_i, va='center', ha='left')
                tile += 1

def plot_convergence_diag(vars):
    """ plot autocorrelation for all stochs in a dict or dict of dicts"""
    # count number of stochastics in model
    cells = 0
    stochs = []
    for k in vars.keys():
        # handle dicts and dicts of dicts by making a list of nodes
        if isinstance(vars[k], dict):
            nodes = vars[k].values()
        else:
            nodes = [vars[k]]

        for n in nodes:
            if isinstance(n, mc.Stochastic) and not n.observed:
                stochs.append(n)
                cells += len(pl.atleast_1d(n.value))

    # for each stoch, make an autocorrelation plot for each dimension
    rows = pl.floor(pl.sqrt(cells))
    cols = pl.ceil(cells/rows)

    tile = 1
    for s in sorted(stochs):
        trace = s.trace()
        if len(trace.shape) == 1:
            trace = trace.reshape((len(trace), 1))
        for d in range(len(pl.atleast_1d(s.value))):
            pl.subplot(rows, cols, tile)
            pl.acorr(pl.atleast_2d(trace)[:, d], normed=True, detrend=pl.mlab.detrend_mean, maxlags=50)
            pl.xticks([])
            pl.yticks([])
            l,r,b,t = pl.axis()
            pl.axis([-.5, r, -.1, 1.1])
            pl.title('\n\n%s[%d]'%(s.__name__, d), va='top', ha='center', fontsize=8)

            tile += 1
    pl.subplots_adjust(0,0,1,1,0,0)
    
    
