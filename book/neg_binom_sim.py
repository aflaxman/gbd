""" Confirm that negative binomial parameter posterior has coverage intended"""


### @export 'setup'
import sys
sys.path += ['..']

import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

iter = 10000
burn = 5000
thin = 1
results = {}
xmax = .07

replicates = 1000

residuals = [[], []]
coverage = [[], []]

### @export 'neg-binom-sim-study'
pi_true = .025
delta_true = 5.
n_pred = 1.e9

for i in range(replicates):
    print '\nsimulation replicate %d' % i
    ## generate simulated data
    n = pl.array(pl.exp(mc.rnormal(10, 1**-2, size=16)), dtype=int)
    k = pl.array(mc.rnegative_binomial(n*pi_true, delta_true), dtype=float)
    r = k/n


    ## setup negative binomial model
    pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
    delta = mc.Uninformative('delta', value=100.)

    @mc.potential
    def obs(pi=pi, delta=delta):
        return mc.negative_binomial_like(r*n, pi*n, delta)

    @mc.deterministic
    def pred(pi=pi, delta=delta):
        return mc.rnegative_binomial(pi*n_pred, delta) / float(n_pred)

    ## fit model
    mc.MCMC([pi, delta, obs, pred]).sample(iter, burn, thin)


    ## record results
    for i, stoch in enumerate([pi, pred]):
        median = stoch.stats()['quantiles'][50]
        residuals[i].append(pi_true - median)
        lb, ub = stoch.stats()['95% HPD interval']
        coverage[i].append(lb <= pi_true <= ub)

### @export 'summarize-results'
bias = {}
rmse = {}
percent_coverage = {}

for i, s in enumerate(['pi', 'pred']):
    bias[s] = pl.mean(residuals[i])
    rmse[s] = pl.rms_flat(residuals[i])
    percent_coverage[s] = pl.mean(coverage[i])

print 'bias', bias
print 'rmse', rmse
print 'percent coverage', percent_coverage

book_graphics.save_json('neg_binom_sim.json', vars())
