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

def plot_prediction_over_time(country, data, predicted, age=0, pred_stats=None, cmap=pl.cm.YlGnBu, alpha=.5, more_data=None):
    """ Plot the predicted values for a specific country as a function of time

    Parameters
    ----------
    country : str, iso3 code
    data : data rec
    predicted : pymc trace
    """
    if pred_stats == None:
        pred_stats = predicted.stats(batches=1)
    i_c = [i for i in range(len(data)) if data.country[i] == country and data.age[i] == age]
    T = len(i_c)
    n = len(predicted.trace())

    # possibly need to sort i_c by data.time[i_c]
    
    # plot jittered trace, to illustrate distribution
    pl.plot(data.year[i_c] + .1*pl.randn(n, T),
            predicted.trace()[:, i_c],
            color=cmap(age/100.), linestyle='', marker='.', alpha=.125, zorder=-1)

    # plot estimated trend
    pl.plot(data.year[i_c], pred_stats['mean'][i_c], zorder=2,
            color=cmap(age/100.), linewidth=3, alpha=alpha)

    # plot 95% HPD
    pl.plot(data.year[i_c], pred_stats['95% HPD interval'][i_c], zorder=2,
            color=cmap(age/100.), linewidth=1, alpha=alpha)

    # overlay data
    dark_color = cmap(age/100.)
    dark_color = (dark_color[0]*.5, dark_color[1]*.5, dark_color[2]*.5)

    pl.plot(data.year[i_c], data.y[i_c], zorder=0,
            linestyle='', marker='o', mew=3, mec=dark_color, color='none', ms=8)
    if more_data != None:
        pl.plot(data.year[i_c], more_data.y[i_c], zorder=1,
                linestyle='', marker='x', mew=3, mec=dark_color, ms=8)

def plot_all_predictions_over_time_for_age(data, predicted, age=0, cmap=pl.cm.YlGnBu, alpha=.5, more_data=None):
    """ Plot the predicted values for a specific country as a function of time

    Parameters
    ----------
    data : data rec
    predicted : pymc trace
    """
    # memorize stats to speed country-specific plots
    pred_stats = predicted.stats(batches=1)

    x_min = min(data.year)
    x_max = max(data.year)
    y_min = round(min(pred_stats['mean']), -1)
    y_max = round(max(pred_stats['mean']), -1)
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
            plot_prediction_over_time(country, data, predicted, age=age,
                                      pred_stats=pred_stats, cmap=cmap, alpha=alpha, more_data=more_data)
            pl.title('\n\n'+country, va='top')

            # set the axis
            pl.axis([x_min-2, x_max+2, 1.2*y_min, 1.2*y_max])
            if jj > 0:
                pl.yticks([])
            else:
                pl.yticks([y_min, 0, y_max])
            if ii < len(regions)-1:
                pl.xticks([])
            else:
                pl.xticks(range(x_min, x_max, 5))

    # set the border width correctly
    pl.subplots_adjust(left=.05, right=.99, bottom=.05, top=.99, wspace=0, hspace=0)

def plot_all_predictions_over_time(data, predicted, cmap=pl.cm.YlGnBu, alpha=.5, more_data=None):
    """ Plot the predicted values for a specific country as a function of time for each age

    Parameters
    ----------
    data : data rec
    predicted : pymc trace
    additional optional parameters, to be described
    """
    for a in pl.unique(data.age):
        print 'plotting for age %s' % a
        plot_all_predictions_over_time_for_age(data, predicted, cmap=cmap, alpha=alpha, more_data=more_data, age=a)

def test(data, mc_dict):
    """ Test plots for all mcmc results in the mc_dict"""

    for mod, mod_mc in mc_dict.items():
        print 'testing plots of mod %s' % mod
        plot_prediction_over_time('USA', data, mod_mc.predicted)  # test single plot of model predictions
        plot_all_predictions_over_time(data, mod_mc.predicted)  # test single plot of model predictions

