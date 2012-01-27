import pylab as pl

import dismod3
dpi=120
quarter_page_params = dict(figsize=(8.5,3), dpi=dpi)
half_page_params = dict(figsize=(8.5, 5.5), dpi=dpi)

width=2
marker_size=5
def plot_age_patterns(dm, region='north_america_high_income', year='2005', sex='male',
                      xticks=[0,25,50,75,100], rate_types='excess-mortality remission incidence prevalence'.split(),
                      yticks=None):
    ages = dm.get_estimate_age_mesh()
    pl.figure(**quarter_page_params)
    pl.subplots_adjust(.1, .175, .98, .875, .275)

    key = dismod3.utils.gbd_key_for('%s', region, year, sex)

    for i, rate_type in enumerate(rate_types):
        pl.subplot(1,len(rate_types),i+1)

        dismod3.plotting.plot_intervals(dm, [d for d in dm.data if dm.relevant_to(d, rate_type, region, year, sex)],
                                        print_sample_size=False, plot_error_bars=False, color='k', alpha=1)

        if rate_type.startswith('prevalence'):
            linestyle='-'
        else:
            linestyle='steps-post-'

        plot_rate(dm, key%rate_type, linestyle=linestyle)

        pl.title('%s' % rate_type.capitalize(), va='bottom')
        l,r,b,t=pl.axis()
        pl.xlabel('Age (Years)')
        pl.xticks(xticks[:-1])
        l,r = xticks[0]-2, xticks[-1]+2

        if isinstance(yticks, dict):
            pl.yticks(yticks[rate_type])
            b,t = yticks[rate_type][0], yticks[rate_type][-1]
            h = t-b
            b -= .05*h
            t += .05*h
        elif isinstance(yticks, list):
            # use the same yticks for each subplot, which means they can be closer together
            if i == 0:
                pl.yticks(yticks)
                pl.ylabel('Rate (Per PY)')
            else:
                pl.yticks(yticks, ['' for y in yticks])
            b,t = yticks[0], yticks[-1]
            h = t-b
            b -= .05*h
            t += .05*h
            pl.subplots_adjust(wspace=.0001)
        pl.axis([l, r, b, t])
            

def plot_rate(dm, key, linestyle='steps-post-'):
    if not isinstance(dm.vars[key]['rate_stoch'].trace, bool):
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
        r = dm.vars[key]['rate_stoch'].value
        pl.plot(dm.get_estimate_age_mesh(), r,
                'w', linewidth=3, linestyle=linestyle)
        pl.plot(dm.get_estimate_age_mesh(), r,
                'k', linewidth=1, linestyle=linestyle)
        pl.plot(dm.get_param_age_mesh()[:-1],
                r[dm.get_param_age_mesh()[:-1]],
                'ko', ms=marker_size, mec='white', zorder=2)


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
                fig_params=half_page_params, subplot_params=dict(bottom=.1, right=.99, top=.95, left=.3), **params):
    sorted_indices = (-r).argsort().argsort()

    se = 1.96*pl.sqrt(r*(1-r)/n)
    ms = pl.minimum(25, pl.sqrt(n) / 10.)

    pl.figure(**fig_params)

    for i in range(len(r)):
        pl.errorbar(r[i], sorted_indices[i]*.5, xerr=[[r[i]-max(0, r[i]-se[i])], [se[i]]],
                    fmt='ks', mew=1, mec='white',
                    ms=5) #ms[i])
        if data_labels:
            pl.text(-2*xmax/50, sorted_indices[i]*.5, data_labels[i], ha='right', va='center')
    pl.yticks([])
    if not data_labels:
        pl.text(-2*xmax/50, (len(sorted_indices)-1)*.25, 'Simulated Study Data', rotation=90, ha='right', va='center')
    if not model_keys:
        if results:
            model_keys = results.keys()
        else:
            model_keys = []

    for i, k in enumerate(model_keys):
        pl.text(-2*xmax/50, -(i+2), k, ha='right', va='center')

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

        pl.errorbar(pi_med, -(i+2), xerr=xerr,
                    fmt='ko', mew=1, mec='white', ms=5)

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

        pl.errorbar(pi_med, -(i+2)+.25, xerr=xerr,
                    fmt='k^', mew=1, mec='white', ms=8)

        pl.hlines([-1], -1, 1, linewidth=1, linestyle='solid', color='k')
        pl.text(-2*xmax/50, -1., 'Model Estimate of Pop. Rate:', va='center', ha='right')


    l,r,b,t=pl.axis()
    b -= .5
    t += .5

    if pi_true:
        pl.vlines([pi_true], b, t, linewidth=1, linestyle='dashed', color='k')
        pl.text(pi_true, t, '\n $\\pi_{true}$', ha='left', va='top')

    pl.axis([-xmax/50., xmax, b, t])
    pl.subplots_adjust(**subplot_params)
    pl.xlabel('Rate (per PY)')

    if fname:
        pl.savefig(fname)
