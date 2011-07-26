import pylab as pl

def plot_age_patterns(dm):
    ages = dm.get_estimate_age_mesh()
    for i in dm.get_param_age_mesh():
        if i > 0:
            ages[i-1] += .5
            ages[i] -= .5
    pl.figure(figsize=(24./1.25, 6./1.25))
    pl.subplots_adjust(.05, .1, .98, .9, .1)

    for i, rate_type in enumerate('prevalence remission incidence excess-mortality'.split()):
        pl.subplot(1,4,4-i)
        pl.plot(ages,
                dm.vars['%s+north_america_high_income+2005+male'%rate_type]['rate_stoch'].value,
                'b', linewidth=2, zorder=1,
                label=rate_type.replace('-', ' ').capitalize())
        if rate_type != 'prevalence':
            # leave off the last point of the param age mesh, since it actually doesn't get used, and hence can be confusing
            pl.plot(dm.get_param_age_mesh()[:-1],
                    pl.exp(dm.vars['%s+north_america_high_income+2005+male'%rate_type]['age_coeffs_mesh'].value)[:-1],
                    'bo', mec='white', zorder=2)

        pl.xlabel('Age (Years)')
        pl.title('%s (Per PY)' % rate_type.capitalize())
        pl.axis([-2, 102, -.01, .39])
        pl.yticks([0.,.1,.2,.3])

        if rate_type == 'excess-mortality':
            pl.semilogy(ages,
                    dm.vars['%s+north_america_high_income+2005+male'%'mortality']['rate_stoch'].value,
                    'r--', linewidth=2,
                    label='With-condition mortality')
            pl.plot(ages,
                    dm.vars['%s+north_america_high_income+2005+male'%'m'].value,
                    'c-.', linewidth=2,
                    label='Background mortality')
            pl.plot(dm.get_estimate_age_mesh(),
                    dm.vars['%s+north_america_high_income+2005+male'%'all-cause_mortality'],
                    'k:', linewidth=2,
                    label='All-cause mortality')
            pl.title('Mortality (Per PY)')
            pl.legend(loc='upper center', ncol=1, fancybox=True)
            pl.axis([-2, 102, 5.e-5, 500.])
