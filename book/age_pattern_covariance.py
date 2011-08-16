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

scale = dict(Very=1000, Moderately=100, Slightly=50)
linestyle = dict(Slightly='solid', Very='dashed', Moderately='dotted')
for col, smoothness in enumerate(['Slightly', 'Moderately', 'Very']):
    mesh = pl.arange(-101, 101)
    C = mc.gp.FullRankCovariance(mc.gp.matern.euclidean,
                                 amp=1.,
                                 scale=scale[smoothness],
                                 diff_degree=2)

    pl.plot(mesh, C(mesh, [0]),
            marker='', color='black', linewidth=1,
            linestyle=linestyle[smoothness],
            zorder=1,
            label='%s ($\\rho = %d$)' % (smoothness, scale[smoothness]))
pl.xticks([-25, 0, 25,50,75])
pl.xlabel('$\\Delta$ Age (Years)')
pl.yticks([0, .25, .5, .75, 1.0])
pl.ylabel('Autocovariance')

pl.axis([-30, 90, -.1, 1.1])

pl.legend(loc='lower left', fancybox=True, shadow=True)

pl.subplots_adjust(left=.1, bottom=.2, top=.95, right=.95, wspace=0)
pl.savefig('smoothness_covariance.png')
pl.savefig('smoothness_covariance.pdf')

pl.show()

### @export 'store-results'

book_graphics.save_json('age_patterns_covariance.json', vars())

