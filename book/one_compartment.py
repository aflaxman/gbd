# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pylab as pl
import scipy.integrate

# <markdowncell>

# One Compartment Model
# =====================
# 
# Constant Flows
# --------------
# 
# When $b$ and $m$ are constant with respect to time (and the model is for a single age group),
# there is a simple closed-form solution to the compartmental model.
# 
# The model is represented by:
# 
#       b  +-----+  m
#     ---->|  S  |---->
#          +-----+
# 
# In equations this represents
# $$S' = b-m,$$
# $$b = h_b S,$$
# $$m = h_m S.$$
# 
# And then if $h_b$ and $h_m$ are constants, the closed form solution is
# $$S = S_0 e^{(h_b-h_m)t}.$$

# <codecell>

t = 1.*pl.arange(101)
S_0 = 1.

h_b = .02*pl.ones_like(t)
h_m = .01*pl.ones_like(t)

# <codecell>

def one_compartment_ode(S, t, h_b, h_m):
    # piecewise-constant functions of time implementend as array
    t = int(pl.clip(t, 0, len(h_b)-1))
    return (h_b[t]-h_m[t])*S

# <markdowncell>

# SciPy has two wrappers for the odepack solver, <code>scipy.integrate.odeint</code> (function-based) and <code>scipy.integrate.ode</code> (object-based).
# 
# The function-based wrapper has less overhead, so I'll use it here.  The object-based wrapper has a nicer format for some applications.

# <codecell>

S_approx = scipy.integrate.odeint(one_compartment_ode, S_0, t, (h_b,h_m))
S_exact = S_0*pl.exp((h_b[0]-h_m[0])*t)

# <codecell>
pl.figure()
pl.plot(t, 100*(S_approx.reshape(101) - S_exact)/S_exact, 'k')
pl.hlines([0],0,100,color='k',linestyle='--')
pl.ylabel('Relative error (%)')
pl.xlabel('Age (years)')
yt = [-2e-6, 0, 2e-6, 4.e-6]
pl.yticks(yt, ['%.6f'%y for y in yt])
pl.axis([-5,105,-2.5e-6, 5e-6])
pl.title('ODE error for constant rates')

# <codecell>

import sys
sys.path += ['..', '../..']

import book_graphics

# <codecell>

def plot_system(panel):
    pl.figure(**book_graphics.quarter_page_params)
    #pl.subplots_adjust(.1, .175, .98, .875, .275, 0)
    b = .03
    
    pl.axes([.2, (.875-.175)/2 + .175 + b, (.9 - .2 - .075) * 1 / 3 , (.875-.175)/2])
    pl.step(t, h_b, 'k-')
    pl.xticks(pl.arange(0, 100, 25), fontsize='x-large')
    pl.yticks([0., .01, .02], fontsize='x-large')
    pl.axis([-5, 105, 0, .03])
    pl.ylabel('$h_b(t)$', rotation='horizontal', fontsize='xx-large')
    
    pl.axes([.2, (.875-.175)*0/2 + .175 + b, (.9 - .2 - .075) * 1 / 3 , (.875-.175)/2])
    pl.step(t, h_m, 'k-')
    pl.yticks([0., .01, .02], fontsize='x-large')
    pl.xticks(pl.arange(0, 100, 25), fontsize='x-large')
    pl.axis([-5, 105, 0, .03])
    pl.xlabel('Age (years)', fontsize='xx-large')
    pl.ylabel('$h_m(t)$', rotation='horizontal', fontsize='xx-large')
    
    pl.axes([.2 + .075 + (.9 - .1) / 3, .175 + b, (.9 - .2 - .075) * 2 / 3 , .875-.175])
    pl.plot(t, S_approx, 'k-')
    pl.xticks(pl.arange(0, 101, 25), fontsize='x-large')
    pl.yticks(pl.arange(.5, 3.5, .5), fontsize='x-large')
    pl.axis([-5, 105, .9, 3.1])
    pl.xlabel('Age (years)', fontsize='xx-large')
    pl.ylabel('$S(t)$', rotation='horizontal', fontsize='xx-large')
    
    if panel == 'b':
        pl.axis([-5,105,.9,1.6])

    pl.axes([0,0,1,1], frameon=False)
    pl.xticks([])
    pl.yticks([])
    pl.figtext(0, 1, '\n (%s)'%panel, ha='left', va='top', fontsize='xx-large')

# <codecell>

plot_system('a')

pl.savefig('one_compartment_constant_rate.pdf')

# <markdowncell>

# Time-varying Rates
# ------------------
# 
# When $b$ and $m$ are not constant with respect to time, the ODE does not have a closed form solution.
# If $b$ and $m$ are piecewise-constant, there an exact solution can be computed iteratively.

# <codecell>

h_b = .01*(2. - t/t.max())
h_m = .005*(2. + t/t.max())

# <codecell>

S_approx = scipy.integrate.odeint(one_compartment_ode, S_0, t, (h_b,h_m))

S_exact = 1. * pl.zeros_like(t)
S_exact[0] = S_0
for i in range(len(t)-1):
    S_exact[i+1] = S_exact[i]*pl.exp((h_b[i]-h_m[i])*(t[i+1]-t[i]))

# <codecell>
pl.figure()
pl.plot(t, 100*(S_approx.reshape(101) - S_exact)/S_exact, 'k')
pl.hlines([0],0,100,color='k',linestyle='--')
pl.ylabel('Relative error (%)')
pl.xlabel('Age (years)')
yt = [-10e-5, 0, 10e-5, 20e-5, 30e-5]
pl.yticks(yt, ['%.4f'%y for y in yt])
pl.axis([-5,105,-10e-5, 30e-5])
pl.title('ODE error for piecewise-constant rates')

# <codecell>

plot_system('b')
pl.savefig('one_compartment_varying_rate.pdf')

# <codecell>

pl.show()
