""" Functions for visualizing the results of fitting the models in model.py

Notes
-----

what is it that we want to explore?  

* posterior predicted values are interesting, and can be compared well
  with the actual values, when available.

* the autocorrelation of the traces are interesting, and function as
  guides to mixing; it might be good to order them by magnitude of lag
  1 correlation, to focus on the worst-mixed

* the time trends by country/region are interesting, since this is the
  actual output that we really care about

* the effect coefficients are intersting, since the economists among
  us think these values are meaningful

this adds up to a lot of data.  it could be a good test case for
mspivot, as a way to explore all data visually and creatively.

attributes:  
* for a row of the input data: time, region, country, autocorrelation, effect coefficient posterior
* for a country: plot over time, effect coefficient posterior
* for an effect: plot over countries and regions
* for a region: plot over times and countries
"""

import pylab as pl

def plot_prediction_over_time(country, data, predicted, color='blue', alpha=.5, more_data=None):
    """ Plot the predicted values for a specific country as a function of time

    Parameters
    ----------
    country : str, iso3 code
    data : data rec
    predicted : pymc trace
    """
    pred_stats = predicted.stats(batches=1)
    i_c = [i for i in range(len(data)) if data.country[i] == country]
    T = len(i_c)
    n = len(predicted.trace())

    # possibly need to sort i_c by data.time[i_c]
    
    # plot jittered trace, to illustrate distribution
    pl.plot(data.year[i_c] + .1*pl.randn(n, T),
            predicted.trace()[:, i_c],
            color='k', linestyle='', marker='.', alpha=.125, zorder=-1)

    # plot estimated trend
    pl.plot(data.year[i_c], pred_stats['mean'][i_c],
            color=color, linewidth=3, alpha=alpha)

    # plot 95% HPD
    pl.plot(data.year[i_c], pred_stats['95% HPD interval'][i_c],
            color=color, linestyle='dashed', linewidth=1, alpha=alpha)

    # overlay data
    if more_data != None:
        pl.plot(data.year[i_c], more_data.y[i_c],
                linestyle='', marker='x', mew=3, mec='g', ms=8)
    pl.plot(data.year[i_c], data.y[i_c],
            linestyle='', marker='x', mew=3, mec='r', ms=8)

def plot_all_predictions_over_time(data, predicted, color='blue', alpha=.5, more_data=None):
    """ Plot the predicted values for a specific country as a function of time

    Parameters
    ----------
    data : data rec
    predicted : pymc trace
    """
    # memorize stats to speed country-specific plots (HACK)
    stats_func = predicted.stats  # save stats function
    stats_val = stats_func(batches=1)  # save results of calling stats function
    predicted.stats = lambda batches=1: stats_val  # replace function that does calculation with function that returns memorized results

    y_min = round(min(predicted.stats()['95% HPD interval'][:,0]), -1)
    y_max = round(max(predicted.stats()['95% HPD interval'][:,1]), -1)
    max_countries = 4   # FIXME: don't hard-code constant 4
    regions = sorted(set(data.region))[:4] # FIXME: don't hard-code constant 4
    for ii, region in enumerate(regions):
        print region

        # label the row
        pl.subplot(len(regions), max_countries, ii*max_countries + 1)
        pl.ylabel(region)

        countries = [data.country[i] for i in range(len(data)) if data.region[i] == region]
        countries = sorted(set(countries))[:max_countries]
        for jj, country in enumerate(countries):
            # plot and label the cell
            pl.subplot(len(regions), max_countries, ii*max_countries + jj + 1)
            plot_prediction_over_time(country, data, predicted, color=color, alpha=alpha, more_data=more_data)
            pl.title('\n\n'+country, va='top')

            # set the axis
            pl.axis([1988, 2007, 1.2*y_min, 1.2*y_max]) # FIXME: don't hard-code constants
            if jj > 0:
                pl.yticks([])
            else:
                pl.yticks([y_min, 0, y_max])
            if ii < len(regions)-1:
                pl.xticks([])
            else:
                pl.xticks([1990, 1995, 2000, 2005]) # FIXME: don't hard-code constants

    # set the border width correctly
    pl.subplots_adjust(left=.05, right=.99, bottom=.05, top=.99, wspace=0, hspace=0)

    # undo memorization hack
    predicted.stats = stats_func


def test(data, mc_dict):
    """ Test plots for all mcmc results in the mc_dict"""

    for mod, mod_mc in mc_dict.items():
        print 'testing plots of mod %s' % mod
        plot_prediction_over_time('USA', data, mod_mc.predicted)  # test single plot of model predictions
        plot_all_predictions_over_time(data, mod_mc.predicted)  # test single plot of model predictions

