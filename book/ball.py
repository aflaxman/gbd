# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pylab as pl
import pymc as mc

# <markdowncell>

# Uniform points in an $n$-dimensional ball
# =========================================
# 
# This notebook implements and compares samplers in PyMC
# to sample uniformly from an $n$-dimensional ball,
# i.e to sample from the set
# $$
# \mathbf{B}_n = \\{x \in \mathbf{R}^n: \|x\|\leq 1\\}
# $$

# <codecell>

mc.np.random.seed(1234567)

# simple model
n = 2
X = [mc.Uninformative('X_%d'%i, value=0) for i in range(n)]
@mc.potential
def in_ball(X=X):
    if X[0]**2 + X[1]**2 <= 1.:
        return 0
    else:
        return -pl.inf

# <codecell>

class UniformBall(mc.Gibbs):
    def __init__(self, stochastic, others, verbose=None):
        self.others = others
        self.conjugate = True # pymc will include a Metropolis rejection step on top of the proposal if this is false
        mc.Gibbs.__init__(self, stochastic, verbose)
    
    def propose(self):
        x_other = [X_i.value for X_i in self.others]
        max_val = pl.sqrt(1. - pl.dot(x_other, x_other))
        self.stochastic.value = mc.runiform(-max_val, max_val)

# <codecell>

m = mc.MCMC([X, in_ball])
for i in range(n):
    m.use_step_method(UniformBall, X[i], [X[j] for j in range(n) if j != i])

# <codecell>

m.sample(100, progress_bar=False)

# <codecell>

def plot_trace(X, scale=1., angle=0.):
    fig = pl.figure(figsize=(12,4.75))
    
    ax1 = fig.add_subplot(1, 2, 1)
    # plot boundary
    t = pl.arange(0,2*pl.pi,.01)
    ax1.plot(pl.cos(angle)*pl.cos(t) - pl.sin(angle)*pl.sin(t)/scale, pl.cos(angle)*pl.sin(t)/scale + pl.sin(angle)*pl.cos(t), 'k:')
    
    # plot samples
    if isinstance(X, mc.Stochastic):
        tr = [X.trace()[:,0], X.trace()[:,1]]
    else:
        tr = [X[0].trace(), X[1].trace()]

    ax1.plot(tr[0], tr[1], 'ko-')
        
    # decorate plot
    pl.xticks(size=18)
    pl.yticks(size=18)
    pl.xlabel('$X_1$', fontsize=32)
    pl.ylabel('$X_2$', fontsize=32, rotation=0)
    pl.axis([-1.1,1.1,-1.1,1.1])
    pl.text(-1,1,'(a)', fontsize=18, va='top', ha='left')

    
    for i in range(2):
        if i == 0:
            ax2 = fig.add_subplot(2, 4, 3+4*i)
            ax2.plot(tr[i], 'k', drawstyle='steps-mid')
        else:
            ax2a = fig.add_subplot(2, 4, 3+4*i, sharex=ax2)
            ax2a.plot(tr[i], 'k', drawstyle='steps-mid')
            pl.xlabel('Sample', fontsize=24)
        
        pl.xticks([25,50,75], size=18)
        pl.yticks([-.5,0,.5], size=18)
        pl.ylabel('$X_%d$'%(i+1), fontsize=32, rotation=0)
        pl.axis([-5,105,-1.5,1.5])
        pl.text(-1,1.25,'(%s)'%'bc'[i], fontsize=18, va='top', ha='left')
        
        if i == 0:
            ax3 = fig.add_subplot(2, 4, 4+4*i)
            ax3.acorr(tr[i].reshape(100), color='k')
        else:
            ax3a = fig.add_subplot(2, 4, 4+4*i, sharex=ax3)
            ax3a.acorr(tr[i].reshape(100), color='k')
            pl.xlabel('Autocorrelation', fontsize=24)
        
        pl.xticks([-5,0,5], size=18)
        pl.yticks([0., .5, 1], size=18)
        pl.axis([-12,12,-.1,1.1])
        pl.text(-10,1,'(%s)'%'de'[i], fontsize=18, va='top', ha='left')
        
    pl.setp(ax2.get_xticklabels(), visible=False)
    pl.setp(ax3.get_xticklabels(), visible=False)
    pl.subplots_adjust(wspace=.55, hspace=.1, bottom=.2,left=.15)

# <codecell>

plot_trace(X, 1, 0.)
pl.savefig('book/graphics/gibbs-ball.pdf')

# <markdowncell>

# Now with the Metropolis sampler
# ---------------------------------

# <codecell>

mc.np.random.seed(123456789)

# <codecell>

# simple model

n = 2
X = mc.Uninformative('X', value=[0,0])
@mc.potential
def in_ball(X=X, s=3., t=pl.pi/4.):
    if (pl.cos(t)*X[0] + pl.sin(t)*X[1])**2 + s**2*(pl.cos(t)*X[1] -pl.sin(t)*X[0])**2 <= 1.:
        return 0
    else:
        return -pl.inf
        
m = mc.MCMC([X, in_ball])

m.sample(100, progress_bar=False)

# <codecell>

plot_trace(X, 3, pl.pi/4)
pl.savefig('book/graphics/metropolis-ball.pdf')

# <markdowncell>

# Now with Adaptive Metropolis

# <codecell>

mc.np.random.seed(1234567)

# simple model
n = 2
X = mc.Uninformative('X', value=[0,0])
@mc.potential
def in_ball(X=X, s=3., t=pl.pi/4):
    if (pl.cos(t)*X[0] + pl.sin(t)*X[1])**2 + s**2*(pl.cos(t)*X[1] -pl.sin(t)*X[0])**2 <= 1.:
        return 0
    else:
        return -pl.inf
        
m = mc.MCMC([X, in_ball])
m.use_step_method(mc.AdaptiveMetropolis, X)

# <codecell>

m.sample(100, progress_bar=False)

plot_trace(X, 3, pl.pi/4)
pl.savefig('book/graphics/am-ball-1.pdf')

# <codecell>

m.sample(iter=20100, burn=20000, progress_bar=False)

plot_trace(X, 3, pl.pi/4)
pl.savefig('book/graphics/am-ball-2.pdf')

pl.show()


