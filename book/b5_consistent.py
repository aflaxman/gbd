""" script to demonstrate how consistent region fits, through inconsistent country fits, yields inconsistent region fits"""

import pandas

X = pandas.read_csv('/home/j/Project/dismod/dismod_status/prod/dm-20084/posterior/dm-20084-prevalence-north_africa_middle_east-male-2005.csv', index_col=None)

Y = pandas.read_csv('/home/j/Project/dismod/dismod_status/prod/dm-19807/posterior/dm-19807-prevalence-north_africa_middle_east-male-2005.csv', index_col=None)

import pylab as pl


def weighted_age(df):                                                    
    return (df.filter(like='Draw').T*df['Population']/df['Population'].sum()).T.sum()


pl.figure()
for iso in list(pl.unique(X['Iso3'])):
    pl.plot(X[X['Iso3']==iso].filter(like='Draw').mean(1).__array__(), label=iso)
pl.semilogy([1],[1])

Z = X.groupby('Age').apply(weighted_age)
plot(Z.mean(1).__array__(), color='red', linewidth=3, alpha=.5, label='Inconsistent NA/ME')

pl.legend()
pl.axis([-5,130,1e-6,2])



pl.figure()
for iso in list(pl.unique(Y['Iso3'])):
    pl.plot(Y[(Y['Iso3']==iso)&(Y['Rate type']=='prevalence')].filter(like='Draw').mean(1).__array__(), label=iso)

pl.semilogy([1],[1])

Z = Y[Y['Rate type'] == 'prevalence'].groupby('Age').apply(weighted_age)
pl.plot(Z.mean(1).__array__(), color='red', linewidth=3, alpha=.5, label='Inconsistent NA/ME')

pl.legend()
pl.axis([-5,130,1e-6,2])




import dismod3
dm = dismod3.load_disease_model(19807)
import fit_posterior
fit_posterior.fit_posterior(dm, 'north_africa_middle_east', 'male', '2005', map_only=True)
X = pandas.read_csv('/var/tmp/dismod_working/test/dm-19807/posterior/dm-19807-north_africa_middle_east-male-2005.csv', index_col=None)
pl.figure()
for iso in list(pl.unique(X['Iso3'])):
    pl.plot(X[(X['Iso3']==iso)].filter(like='Draw').mean(1).__array__(), label=iso)

pl.semilogy([1],[1])


Z = X.groupby('Age').apply(weighted_age)
plot(Z.mean(1).__array__(), color='red', linewidth=3, alpha=.5, label='Inconsistent NA/ME')

plot(dm.vars['prevalence+north_africa_middle_east+2005+male']['rate_stoch'].stats()['mean'], color='red', linewidth=3, alpha=.5, label='Mean of Consistent NA/ME')


pl.legend()
pl.axis([-5,130,1e-6,2])



dm = dismod3.load_disease_model(19807)
import fit_emp_prior
fit_emp_prior.fit_emp_prior(dm, 'prevalence')
fit_posterior.fit_posterior(dm, 'north_africa_middle_east', 'male', '1990')

pl.plot(dm.get_mcmc('mean', 'prevalence+north_africa_middle_east+1990+male'), color='black', linewidth=3, alpha=.5, label='Mean of Inconsistent NA/ME')
pl.plot(dm.vars['prevalence+north_africa_middle_east+1990+male']['rate_stoch'].stats()['mean'], color='red', linewidth=3, alpha=.5, label='Mean of Consistent NA/ME')
pl.semilogy([1],[1])

Y = pandas.read_csv('/var/tmp/dismod_working/test/dm-19807/posterior/dm-19807-prevalence-north_africa_middle_east-male-1990.csv', index_col=None)
for iso in list(pl.unique(Y['Iso3'])):
    pl.plot(Y[Y['Iso3']==iso].filter(like='Draw').mean(1).__array__(), label=iso)

pl.legend()
pl.axis([-5,130,1e-6,2])
