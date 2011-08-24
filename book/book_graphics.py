import pylab as pl

import dismod3
dpi=120
quarter_page_params = dict(figsize=(8.5,3), dpi=dpi)
half_page_params = dict(figsize=(8.5, 5.5), dpi=dpi)

large=20
width=2
marker_size=5
def plot_age_patterns(dm, r='north_america_high_income', y='2005', s='male'):
    ages = dm.get_estimate_age_mesh()
    for i in dm.get_param_age_mesh():
        if i > 0:
            ages[i-1] += .5
            ages[i] -= .5
    pl.figure(figsize=(24./1.25, 6./1.25))
    pl.subplots_adjust(.05, .175, .98, .875, .175)

    key = dismod3.utils.gbd_key_for('%s', r, y, s)

    for i, rate_type in enumerate('prevalence remission incidence excess-mortality'.split()):
        pl.subplot(1,4,4-i)
        if rate_type == 'prevalence':
            pl.plot(range(dismod3.settings.MAX_AGE),
                    dm.vars[key%rate_type]['rate_stoch'].value,
                    'k-', linewidth=width)
            pl.plot(dm.get_param_age_mesh()[:-1],
                    dm.vars[key%rate_type]['rate_stoch'].value[dm.get_param_age_mesh()[:-1]],
                    'ko', ms=marker_size, mec='white', zorder=2)
        else:
            pl.plot(ages,
                    dm.vars[key%rate_type]['rate_stoch'].value,
                    'k', linewidth=width, zorder=1,
                    label=rate_type.replace('-', ' ').capitalize())
            # leave off the last point of the param age mesh, since it actually doesn't get used, and hence can be confusing
            pl.plot(dm.get_param_age_mesh()[:-1],
                    pl.exp(dm.vars[key%rate_type]['age_coeffs_mesh'].value)[:-1],
                    'ko', ms=marker_size, mec='white', zorder=2)
        dismod3.plotting.plot_intervals(dm, [d for d in dm.data if dm.relevant_to(d, rate_type, 'all', 'all', 'all')],
                                        print_sample_size=False, plot_error_bars=False, alpha=1)

        pl.title('%s (Per PY)' % rate_type.capitalize(), fontsize=large)
        pl.axis([-2, 102, -.01, .39])
        pl.yticks([0.,.1,.2,.3], fontsize=large)

        if rate_type == 'excess-mortality':
            pl.semilogy(ages,
                    dm.vars[key%'mortality']['rate_stoch'].value,
                    'r--', linewidth=width,
                    label='With-condition mortality')
            pl.plot(ages,
                    dm.vars[key%'m_background'].value,
                    'c-.', linewidth=width,
                    label='Background mortality')
            pl.plot(dm.get_estimate_age_mesh(),
                    dm.vars[key%'all-cause_mortality'],
                    'k:', linewidth=width,
                    label='All-cause mortality')
            pl.title('Mortality (Per PY)', fontsize=large)
            pl.legend(loc='upper center', ncol=1, fancybox=True)
            #pl.axis([-2, 102, 5.e-5, 500.])
            pl.yticks(fontsize=large)
        pl.xlabel('Age (Years)', fontsize=large)
        pl.xticks([0,25,50,75],fontsize=large)



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
        pl.errorbar(r[i], sorted_indices[i]*.5, xerr=se[i],
                    fmt='gs', mew=1, mec='white',
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
                    fmt='bo', mew=1, mec='white', ms=5)

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
                    fmt='r^', mew=1, mec='white', ms=5)

        pl.hlines([-1], -1, 1, linewidth=1, linestyle='solid', color='k')
        pl.text(-2*xmax/50, -1., 'Model Estimate of Pop. Rate:', va='center', ha='right')


    l,r,b,t=pl.axis()
    b -= .5
    t += .5

    if pi_true:
        pl.vlines([pi_true], b, t, linewidth=1, linestyle='dashed', color='r')
        pl.text(pi_true, t, '\n $\\pi_{true}$', ha='left', va='top')

    pl.axis([-xmax/50., xmax, b, t])
    pl.subplots_adjust(**subplot_params)
    pl.xlabel('Rate (per PY)')

    if fname:
        pl.savefig(fname)
