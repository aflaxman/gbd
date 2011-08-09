import pylab as pl

import dismod3
dpi=120
quarter_page_params = dict(figsize=(8.5,3), dpi=dpi)
half_page_params = dict(figsize=(8.5, 5.5), dpi=dpi)

large = 24
width=5
marker_size=10
def plot_age_patterns(dm):
    ages = dm.get_estimate_age_mesh()
    for i in dm.get_param_age_mesh():
        if i > 0:
            ages[i-1] += .5
            ages[i] -= .5
    pl.figure(figsize=(24./1.25, 6./1.25))
    pl.subplots_adjust(.05, .175, .98, .875, .175)

    for i, rate_type in enumerate('prevalence remission incidence excess-mortality'.split()):
        pl.subplot(1,4,4-i)
        if rate_type == 'prevalence':
            pl.plot(range(dismod3.settings.MAX_AGE),
                    dm.vars['%s+north_america_high_income+2005+male'%rate_type]['rate_stoch'].value,
                    'b--', linewidth=width)
            pl.plot(dm.get_param_age_mesh()[:-1],
                    dm.vars['%s+north_america_high_income+2005+male'%rate_type]['rate_stoch'].value[dm.get_param_age_mesh()[:-1]],
                    'bo', ms=marker_size, mec='white', zorder=2)
        else:
            pl.plot(ages,
                    dm.vars['%s+north_america_high_income+2005+male'%rate_type]['rate_stoch'].value,
                    'b', linewidth=width, zorder=1,
                    label=rate_type.replace('-', ' ').capitalize())
            # leave off the last point of the param age mesh, since it actually doesn't get used, and hence can be confusing
            pl.plot(dm.get_param_age_mesh()[:-1],
                    pl.exp(dm.vars['%s+north_america_high_income+2005+male'%rate_type]['age_coeffs_mesh'].value)[:-1],
                    'bo', ms=marker_size, mec='white', zorder=2)

        pl.title('%s (Per PY)' % rate_type.capitalize(), fontsize=large)
        pl.axis([-2, 102, -.01, .39])
        pl.yticks([0.,.1,.2,.3], fontsize=large)

        if rate_type == 'excess-mortality':
            pl.semilogy(ages,
                    dm.vars['%s+north_america_high_income+2005+male'%'mortality']['rate_stoch'].value,
                    'r--', linewidth=width,
                    label='With-condition mortality')
            pl.plot(ages,
                    dm.vars['%s+north_america_high_income+2005+male'%'m'].value,
                    'c-.', linewidth=width,
                    label='Background mortality')
            pl.plot(dm.get_estimate_age_mesh(),
                    dm.vars['%s+north_america_high_income+2005+male'%'all-cause_mortality'],
                    'k:', linewidth=width,
                    label='All-cause mortality')
            pl.title('Mortality (Per PY)', fontsize=large)
            pl.legend(loc='upper center', ncol=1, fancybox=True)
            pl.axis([-2, 102, 5.e-5, 500.])
            pl.yticks(fontsize=large)
        pl.xlabel('Age (Years)', fontsize=large)
        pl.xticks([0,25,50,75],fontsize=large)
