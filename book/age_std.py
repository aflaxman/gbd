""" Compare age-standardizing model to midpoint approx"""


### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)


### @export 'initialize-model'


### @export 'initial-rates'

book_graphics.plot_age_patterns(dm)
pl.savefig('initial.png')
