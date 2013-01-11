import sys
sys.path += ['../gbd', '../gbd/book', '../dm3-computation_only/', '../dm3-computation_only/book']
import pylab as pl
import pymc as mc
import pandas

import dismod3
reload(dismod3)

import book_graphics
reload(book_graphics)
import matplotlib.pyplot as plt

# set font
book_graphics.set_font()

def my_axis(ymax):
    pl.axis([10.1,63,-ymax/10.,ymax])
    
	
# def subtitle(s):
    # """ title where the panel names appear within each panel"""
    # l,r,b,t=pl.axis()
    # x = l + (r-l)*.05
    # y = t - (t-b)*.05
    # pl.text(x, y, s, ha='left', va='top', size=16)

# def subtitle_third(s):
    # """ title where the panel names appear within each panel"""
    # l,r,b,t=pl.axis()
    # x = l + (r-l)*.05
    # y = t - (t-b)*.2
    # pl.text(x, y, s, ha='left', va='top', size=16)

# def subtitle_fourth(s):
    # """ title where the panel names appear within each panel"""
    # l,r,b,t=pl.axis()
    # x = l + (r-l)*.05
    # y = t - (t-b)*.25
    # pl.text(x, y, s, ha='left', va='top', size=16)
    
def load_new_model():
    try:
        orig_model = dismod3.data.load('/home/j/Project/dismod/output/dm-32404') 
    except:
        model = dismod3.data.load('/home/j/Project/dismod/dismod_status/prod/dm-32404')
    orig_model.keep(areas=['europe_western'])
    return orig_model

best_model = load_new_model()
output = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_pms.csv')

# figure pms-data
pl.figure(**book_graphics.full_page_params)

dismod3.graphics.plot_data_bars(best_model.get_data('p'))
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .25, .5, .75, 1], [0, 25, 50, 75, 100])
my_axis(1.05)

pl.savefig('book/graphics/pms-data.pdf')
pl.savefig('book/graphics/pms-data.png')

# figure pms-priors
import matplotlib.pyplot as plt
fig = pl.figure(**book_graphics.full_page_params)
ymax=.35

ax1 = fig.add_subplot(4,1,1)
#dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')
ax1.plot(pl.array(output['p_unr']), 'k-', linewidth=3, label='$p(a)$ unrestricted')
pl.axis([10.1,63,-ymax/10.,ymax])
pl.yticks([0, .10, .20, .3], [0, 10, 20, 30])
pl.text(62, .31, '$p(a)$ unrestricted', ha='right', va='top',size=16)
book_graphics.subtitle_fourth('(a)')

ax2 = fig.add_subplot(4,1,2, sharex=ax1)
#dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')
ax2.plot(pl.array(output['p_u15']), 'k-', linewidth=3, label='$p(a)=$ $0$ for $a <$ $15$')
pl.axis([10.1,63,-ymax/10.,ymax])
pl.yticks([0, .10, .20, .3], [0, 10, 20, 30])
pl.text(62, .31, '$p(a)=$ $0$, for $a <$ $15$', ha='right', va='top',size=16)
book_graphics.subtitle_fourth('(b)')

ax3 = fig.add_subplot(4,1,3, sharex=ax1)
#dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')
ax3.plot(pl.array(output['p_o50']), 'k-', linewidth=3, label='$p(a)=$ $0$, for $a >$ $50$')
pl.axis([10.1,63,-ymax/10.,ymax])
pl.yticks([0, .10, .20, .3], [0, 10, 20, 30])
pl.text(62, .31, '$p(a)=$ $0$, for $a >$ $50$', ha='right', va='top',size=16)
pl.ylabel('                        Prevalence (%)'+'\n\n', ha='center')
book_graphics.subtitle_fourth('(c)')

ax4 = fig.add_subplot(4,1,4, sharex=ax1)
#dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')
ax4.plot(pl.array(output['p_u15o50']), 'k-', linewidth=3, label='$p(a)=$ $0$, for $a<$ $15$ and $a>$ $50$')
pl.axis([10.1,63,-ymax/10.,ymax])
pl.yticks([0, .10, .20, .3], [0, 10, 20, 30])
pl.text(62, .31, '$p(a)=$ $0$, for $a<$ $15$ and $a>$ $50$', ha='right', va='top',size=16)
book_graphics.subtitle_fourth('(d)')

pl.xlabel('Age (years)')

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)

pl.subplots_adjust(hspace=.1)

pl.savefig('book/graphics/pms-priors.pdf')
pl.savefig('book/graphics/pms-priors.png')

# figure pms-knot_location

#subplots_adjust(hspace=.35)
#left  = 0.125  # the left side of the subplots of the figure
#  right = 0.9    # the right side of the subplots of the figure
#  bottom = 0.1   # the bottom of the subplots of the figure
#  top = 0.9      # the top of the subplots of the figure
#  wspace = 0.2   # the amount of width reserved for blank space between subplots
#  hspace = 0.2 
# l t w h

fig = pl.figure(**book_graphics.full_page_params)
ymax=.45

ax1 = fig.add_subplot(3,1,1)
#dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')
ax1.plot(pl.array(output['k1_m']), 'k-', linewidth=3, label='{32}') #label='0, 15, 32, 50, 100')
ax1.plot(pl.array(output['k1_l']), 'k--', linewidth=3, label='{20}') #label='0, 15, 20, 50, 100')
ax1.plot(pl.array(output['k1_r']), 'k:', linewidth=3, label='{45}') #label='0, 15, 45, 50, 100')
pl.yticks([0, .1, .2, .3, .4], [0, 10, 20, 30, 40])
pl.legend(loc='upper right', fancybox=True, shadow=True, title='Knots at 15,50 and:')
pl.axis([10.1,63,-ymax/10.,ymax])
book_graphics.subtitle_third('(a)')

ax2 = fig.add_subplot(3,1,2, sharex=ax1)
#dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')
ax2.plot(pl.array(output['k2_m']), 'k-', linewidth=3, label='{27, 38}') #label='0, 15, 27, 38, 50, 100')
ax2.plot(pl.array(output['k2_o']), 'k--', linewidth=3, label='{20, 45}') #label='0, 15, 20, 45, 50, 100')
ax2.plot(pl.array(output['k2_i']), 'k:', linewidth=3, label='{30, 35}') #label='0, 15, 30, 35, 50, 100')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .1, .2, .3, .4], [0, 10, 20, 30, 40])
pl.legend(loc='upper right', fancybox=True, shadow=True, title='Knots at 15,50 and:')
pl.axis([10.1,63,-ymax/10.,ymax])
book_graphics.subtitle_third('(b)')

ax3 = fig.add_subplot(3,1,3, sharex=ax1)
#dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')
ax3.plot(pl.array(output['k3_m']), 'k-', linewidth=3, label='{23, 32, 41}') #label='0, 15, 23, 32, 41, 50, 100')
ax3.plot(pl.array(output['k3_o']), 'k--', linewidth=3, label='{18, 32, 47}') #label='0, 15, 18, 32, 47, 50, 100')
ax3.plot(pl.array(output['k3_i']), 'k:', linewidth=3, label='{29, 32, 35}') #label='0, 15, 29, 32, 35, 50, 100')
pl.yticks([0, .1, .2, .3, .4], [0, 10, 20, 30, 40])
pl.legend(loc='upper right', fancybox=True, shadow=True, title='Knots at 15,50 and:')
pl.axis([10.1,63,-ymax/10.,ymax])
book_graphics.subtitle_third('(c)')

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
pl.xlabel('Age (years)')
pl.subplots_adjust(hspace=.1)

pl.savefig('book/graphics/pms-knot_location.pdf')
pl.savefig('book/graphics/pms-knot_location.png')

# figure pms-direction
pl.figure(**book_graphics.full_page_params)
dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')
pl.plot(pl.array(output['decreasing']), 'k:', linewidth=3, label='Decreasing prior')
pl.plot(pl.array(output['increasing']), 'k--', linewidth=3, label='Increasing prior')
pl.plot(pl.array(output['priorless']), 'k-', linewidth=3, label='No prior')
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .1, .2, .3, .4], [0, 10, 20, 30, 40])
my_axis(.45)
pl.legend(loc='upper right', fancybox=True, shadow=True)

pl.savefig('book/graphics/pms-direction.pdf')
pl.savefig('book/graphics/pms-direction.png')

# figure pms-best_model
pl.figure(**book_graphics.full_page_params)
x = best_model.parameters['p']['parameter_age_mesh']

dismod3.graphics.plot_data_bars(best_model.get_data('p'), color='grey')
pl.plot(pl.array(output['suggestion']), 'k-', linewidth=3, label='Posterior mean')
pl.plot(x, pl.array(output['suggestion_l'])[x], 'k-', linewidth=1, label='95% HPD interval')
pl.plot(x, pl.array(output['suggestion_u'])[x], 'k-', linewidth=1)
pl.xlabel('Age (years)')
pl.ylabel('Prevalence (%)')
pl.yticks([0, .25, .5, .75, 1], [0, 25, 50, 75, 100])
my_axis(1.05)
pl.legend(loc='upper right', fancybox=True, shadow=True)

pl.savefig('book/graphics/pms-best_model.pdf')
pl.savefig('book/graphics/pms-best_model.png')

pl.show()