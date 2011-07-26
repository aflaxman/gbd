import pylab as pl

def plot_age_patterns(dm):
    ages = dm.get_estimate_age_mesh()
    for i in dm.get_param_age_mesh():
        if i > 0:
            ages[i-1] += .5
            ages[i] -= .5
    pl.figure(figsize=(24,6))
    pl.subplots_adjust(.05, .1, .95, .98, .25)

    for i, rate_type in enumerate('incidence remission excess-mortality prevalence'.split()):
        pl.subplot(1,4,i+1)
        pl.plot(ages,
                dm.vars['%s+north_america_high_income+2005+male'%rate_type]['rate_stoch'].value,
                'b', linewidth=2,
                label=rate_type.replace('-', ' ').capitalize())
        if rate_type != 'prevalence':
            pl.plot(dm.get_param_age_mesh(),
                    pl.exp(dm.vars['%s+north_america_high_income+2005+male'%rate_type]['age_coeffs_mesh'].value),
                    'bo', mec='white')
        pl.xlabel('Age (Years)')
        pl.ylabel('%s (Per PY)' % rate_type.capitalize())

        pl.axis([-2,102,-0.01,.41])

        if rate_type == 'excess-mortality':
            pl.semilogy(ages,
                    dm.vars['%s+north_america_high_income+2005+male'%'mortality']['rate_stoch'].value,
                    'r--', linewidth=2,
                    label='With-condition mortality')
            pl.semilogy(ages,
                    dm.vars['%s+north_america_high_income+2005+male'%'m'].value,
                    'c-.', linewidth=2,
                    label='Background mortality')
            pl.plot(dm.get_estimate_age_mesh(),
                    dm.vars['%s+north_america_high_income+2005+male'%'all-cause_mortality'],
                    'k:', linewidth=2,
                    label='All-cause mortality')
            pl.ylabel('Mortality (Per PY)')
            pl.legend(loc='lower right')
            pl.axis([-2, 102, 1e-4, .41])
