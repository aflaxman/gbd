""" Beginnning of a project to visualize pandas dataframes with the d3
javascript library"""

import pylab as pl
import pandas
import random
import colorbrewer
colors = ['#%x%x%x' % col for col in colorbrewer.Set1[6]]

def scatter(df, x, y, size=None, color=None, groupby=None):
    """ Generate scatter plot

    Parameters
    ----------
    df : pandas.DataFrame
    x, y : str
      columns to use for x- and y-axis
    size, color : str, optional
      columns to use for size and color of scatter
    groupby : str or list, optional
      column or columns to group plots by, generating subplots for
      each member of the grouping

    Results
    -------

    (Eventually) Return a str full of html/javascript that shows this
    scatter when loaded into a web browser

    For now, just use matplotlib
    """
    for col in [x, y, size, color]:
        if col != None:
            assert col in df, 'Column "%s" not found in DataFrame' % col  # TODO: say which param has the bad col name
    # TODO: check that groupby appears and there are not too many groups

    if groupby != None:
        groups = df.groupby(groupby)

        n = len(groups)
        c = pl.ceil(pl.sqrt(n))
        r = pl.ceil(n/c)

        prev_subplot = None
        for i, (g, sub_df) in enumerate(groups):
            prev_subplot = pl.subplot(r, c, i+1, sharex=prev_subplot, sharey=prev_subplot)
            pl.title('\n\n%s = %s' % (groupby, g), fontsize='small', verticalalignment='top')
            scatter(sub_df, x, y, size, color)

            if i == (r-1)*c:
                pl.xlabel(x, ha='left')
            else:
                pl.xlabel('')

            if i == 0:
                pl.ylabel(y, va='top')
            else:
                pl.ylabel('')

        pl.yticks([])
        pl.xticks([])
        pl.subplots_adjust(wspace=0, hspace=0)

    else:
        if size == None:
            s = 100
        else:
            s = 10 + 500*(df[size] - df[size].min()) / (df[size].max() - df[size].min())
            s[pl.isnan(s)] = 100

        if color == None or len(df) == 1:
            c = colors[0]
        else:
            cats = pl.unique(df[color])
            num_cats = float(len(cats))
            cats = pandas.Series(pl.arange(num_cats)/num_cats, index=cats)
            c = df[color].map(cats)

            # TODO: make legend showing color<->category mapping
            # Requirements
            #   Show category name and color
            #   Show marker size and number
            #   Mouse-over to highlight all markers of that color, or near that size
            #   Select only certain parts of the data
            #   Select marker to see all data associated with it

        pl.scatter(jitter(df, x).__array__(), jitter(df, y).__array__(), s=s, c=c, linewidths=0, alpha=.5)
        pl.xlabel(x)
        pl.ylabel(y)

def jitter(df, x, pct=.01):
    """ Jitter column x by a certain percent
    Results
    -------
    Return a pandas.Series with jittered values
    """

    return df[x] + pl.randn(len(df.index)) * pct * df[x].mean()
