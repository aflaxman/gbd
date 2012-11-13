""" Explore features of age pattern priors"""


### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

pl.figure(**book_graphics.quarter_page_params)
### @export 'models-of-varying-smoothness'

from dismod3.utils import rho
linestyle = dict(slightly='solid', very='dashed', moderately='dotted')
for col, smoothness in enumerate(['slightly', 'moderately', 'very']):
    mesh = pl.arange(-101, 101, .1)
    C = mc.gp.FullRankCovariance(mc.gp.matern.euclidean,
                                 amp=1.,
                                 scale=rho[smoothness],
                                 diff_degree=2)

    pl.plot(mesh, C(mesh, [0]),
            marker='', color='black', linewidth=1,
            linestyle=linestyle[smoothness],
            zorder=1,
            label='%s ($\\rho = %d$)' % (smoothness.capitalize(), rho[smoothness]))
pl.xticks([-25, 0, 25,50,75])
pl.xlabel('$\\Delta$ Age (Years)')
pl.yticks([0, .25, .5, .75, 1.0])
pl.ylabel('Autocovariance')

pl.axis([-30, 90, -.1, 1.1])

pl.legend(loc='upper right', fancybox=True, shadow=True)

pl.subplots_adjust(left=.1, bottom=.2, top=.95, right=.95, wspace=0)
pl.savefig('smoothness_covariance.png')
pl.savefig('smoothness_covariance.pdf')

### @export 'store-results'

book_graphics.save_json('age_pattern_covariance.json', vars())

