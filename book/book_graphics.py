import pylab as pl
import matplotlib.mpl as mpl

import dismod3
dpi=120
quarter_page_params = dict(figsize=(10,3), dpi=dpi)
half_page_params = dict(figsize=(11, 4.25), dpi=dpi)
three_quarter_page_params = dict(figsize=(11,5), dpi=dpi)
full_page_params = dict(figsize=(11, 8.5), dpi=dpi)

width=2
marker_size=5

def set_font():
# make all fonts bigger, etc
    mpl.rcParams['axes.titlesize'] = 'x-large'
    mpl.rcParams['axes.labelsize'] = 'x-large'
    mpl.rcParams['xtick.labelsize'] = 'x-large'
    mpl.rcParams['ytick.labelsize'] = 'x-large'
    mpl.rcParams['legend.fancybox'] = True
    mpl.rcParams['legend.fontsize'] = 'x-large'
    mpl.rcParams['text.fontsize'] = 12

def plot_age_patterns(model, region='north_america_high_income', year='2005', sex='male',
                      xticks=[0,25,50,75,100], types='i r f p'.split(),
                      yticks=None,
                      panel=None):
    ages = model.parameters['ages']
    pl.figure(**quarter_page_params)
    pl.subplots_adjust(.1, .175, .98, .98, .5)

    for i, rate_type in enumerate(types):
        if types == 'i r m f p'.split():
            if rate_type == 'p':
                pl.subplot(1, 2, 2)
            else:
                pl.subplot(2, 4, 1 + (i%2) + 4*(i/2))
        else:
            pl.subplot(1, len(types), i+1)
        
        #pl.subplots_adjust(wspace=1.1)

        if rate_type == 'p':
            rate_name = '$\\frac{C}{S+C}$'
            xticks=[0,25,50,75,100]
        else:
            rate_name = '$h_%s$'%rate_type

        plot_rate(model.vars[rate_type])

        pl.xlabel('Age (years)', fontsize='xx-large')

        l,r,b,t=pl.axis()
        if isinstance(yticks, list): pl.xticks(xticks[:-1], fontsize='x-large')
        else: pl.xticks(xticks, fontsize='x-large')
        l,r = xticks[0]-2, xticks[-1]+2

        if isinstance(yticks, dict):
            pl.yticks(yticks[rate_type], fontsize='x-large')
            if rate_type in 'ir':
                pl.xticks(xticks[:-1], ['' for _ in xticks[:-1]], fontsize='x-large')
                pl.xlabel('')
            b,t = yticks[rate_type][0], yticks[rate_type][-1]
            h = t-b
            b -= .05*h
            t += .15*h
            pl.text(l,t,'\n %s' % rate_name, ha='left', va='top', rotation='horizontal', fontsize='xx-large')
            pl.subplots_adjust(bottom=.2,hspace=.1, wspace=.3)
        elif isinstance(yticks, list):
            # use the same yticks for each subplot, which means they can be closer together
            if i == 0:
                pl.yticks(yticks, fontsize='x-large')
                pl.ylabel('Rate (per 1)', fontsize='xx-large')
            else:
                pl.yticks(yticks, ['' for y in yticks])
            b,t = yticks[0], yticks[-1]
            h = t-b
            b -= .05*h
            t += .05*h
            pl.subplots_adjust(bottom=.2,wspace=.0001)
            pl.text(l,t,'\n %s'%rate_name, ha='left', va='top', fontsize='xx-large')
        pl.axis([l, r, b, t])
            
    if panel:
        pl.axes([0,0,1,1],frameon=False)
        pl.xticks([])
        pl.yticks([])
        pl.figtext(0,1,'\n (%s)'%panel, va='top', ha='left', fontsize='xx-large')
        pl.subplots_adjust(left=.11)


def plot_rate(vars):
    if not isinstance(vars['mu_age'].trace, bool):
        for r in dm.vars[key]['rate_stoch'].trace()[::10]:
            pl.plot(dm.get_estimate_age_mesh(), r, '-', color='grey', linewidth=2, zorder=-100, linestyle=linestyle)
        r = dm.vars[key]['rate_stoch'].stats()['quantiles'][50]
        pl.plot(dm.get_estimate_age_mesh(), r,
                linewidth=3, color='white', linestyle=linestyle)
        pl.plot(dm.get_estimate_age_mesh(), r,
                linewidth=1, color='black', linestyle=linestyle)
        pl.plot(dm.get_param_age_mesh()[:-1],
                r[dm.get_param_age_mesh()[:-1]],
                'ko', ms=marker_size, mec='white', zorder=2)
    else:
        r = vars['mu_age'].value
        x = pl.arange(len(r))
        pl.plot(x, r,
                'w', linewidth=3)
        pl.plot(x, r,
                'k', linewidth=1)
        if 'knots' in vars:
            a = vars['knots']
            pl.plot(a, r[a], 'ks', 
                    ms=marker_size, mec='white', zorder=2)

def save_json(fname, vars):
    def array_to_list(x):
        try:
            return list(x)
        except:
            return None

    import simplejson as json
    json.dump(vars, open(fname, 'w'),
              indent=2, skipkeys=True, default=array_to_list)



def forest_plot(r, n, pi_true=None, results=None, model_keys=None, data_labels=None, fname=None, xmax=.05,
                fig_params=half_page_params, subplot_params=dict(bottom=.1, right=.99, top=.95, left=.33), **params):
    sorted_indices = (-r).argsort().argsort()

    se = 1.96*pl.sqrt(r*(1-r)/n)
    ms = pl.minimum(25, pl.sqrt(n) / 10.)
    pl.figure(**fig_params)

    for i in range(len(r)):
        pl.errorbar(r[i], sorted_indices[i]*.5-.25, xerr=[[r[i]-max(0, r[i]-se[i])], [se[i]]],
                    fmt='ks', mew=1, mec='white',
                    ms=5) #ms[i])
        if data_labels:
            pl.text(-2*xmax/50, sorted_indices[i]*.5-.25, data_labels[i], ha='right', va='center', fontsize='x-large')
            pl.text(-2*xmax/50, len(sorted_indices)*.5-.25, 'Input data:', va='center', ha='right', fontsize='x-large')
    pl.yticks([])
    pl.xticks(size='large')
    if not data_labels:
        pl.text(-2*xmax/50, (len(sorted_indices)-1)*.25, 'Simulated Study Data', rotation=90, ha='right', va='center', fontsize='x-large')
    if not model_keys:
        if results:
            model_keys = results.keys()
        else:
            model_keys = []

    for i, k in enumerate(model_keys):
        if k == 'Beta binomial':
            k1 = 'Beta-binomial'
            pl.text(-2*xmax/50, -(i*.5+1.75), k1, ha='right', va='center', fontsize='x-large')
        elif k == 'Negative binomial':
            k1 = 'Negative-binomial'
            pl.text(-2*xmax/50, -(i*.5+1.75), k1, ha='right', va='center', fontsize='x-large')
        else: pl.text(-2*xmax/50, -(i*.5+1.75), k, ha='right', va='center', fontsize='x-large')

        # plot prediction posterior
        if '50' in results[k]['pred']['quantiles']: # number becomes string when read back from disk
            pi_med = results[k]['pred']['quantiles']['50']
        else:
            pi_med = results[k]['pred']['quantiles'][50]
        pi_lb = results[k]['pred']['95% HPD interval'][0]
        pi_ub = results[k]['pred']['95% HPD interval'][1]
        n = pi_med*(1-pi_med) / ((pi_ub - pi_lb)/(2*1.96))**2
        xerr = [
            [pi_med - pi_lb],
            [pi_ub - pi_med]
            ]

        #if i == 0:
        #    label = 'Predicted Study Value'
        #else:
        #    label = '_nolabel_'
        #pl.errorbar(pi_med, -(i+2), xerr=xerr,
        #            fmt='ko', mew=1, mec='white', ms=5, label=label)

        # plot parameter posterior
        if '50' in results[k]['pi']['quantiles']: # number becomes string when read back from disk
            pi_med = results[k]['pi']['quantiles']['50']
        else:
            pi_med = results[k]['pi']['quantiles'][50]
        pi_lb = results[k]['pi']['95% HPD interval'][0]
        pi_ub = results[k]['pi']['95% HPD interval'][1]
        n = pi_med*(1-pi_med) / ((pi_ub - pi_lb)/(2*1.96))**2
        xerr = [
            [pi_med - pi_lb],
            [pi_ub - pi_med]
            ]

        if i == 0:
            label = 'Parameter value'
        else:
            label = '_nolabel_'
        pl.errorbar(pi_med, -(i*.5+2)+.25, xerr=xerr,
                    fmt='k^', mew=1, mec='white', ms=8, label=label)

        pl.hlines([-.75], -1, 1, linewidth=1, linestyle='dotted', color='k', label='_nolegend_')
        pl.text(-2*xmax/50, -1.25, 'Model estimate of pop. rate:', va='center', ha='right', fontsize='x-large')

        #pl.legend(loc='lower right', shadow=True, fancybox=True, numpoints=1)

    l,r,b,t=pl.axis()
    b -= .5
    t += .75

    if pi_true:
        pl.vlines([pi_true], b, t, linewidth=1, linestyle='dashed', color='k')
        pl.text(pi_true, t, '\n $\\pi_{true}$', ha='left', va='top', fontsize='xx-large')

    pl.axis([-xmax/50., xmax, b, t])
    pl.subplots_adjust(**subplot_params)
    pl.xlabel('Rate (per PY)', fontsize='x-large')

    if fname:
        pl.savefig(fname)
