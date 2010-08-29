""" Functions for generating synthetic data"""

def generate_fe():
    """ replace data.csv with random data based on a fixed effects model

    This function generates data for all countries in all regions, based on the model::

        Y_r,c,t = beta . X_r,c,t + e_r,c,t
        e_r,c,t ~ N(0,1)

    """
    import csv
    from pylab import randn, dot

    c4 = dict([[d[0], d[1:]] for d in csv.reader(open('../country_region.csv'))])
    c4.pop('World')

    data = [['region', 'country', 'year', 'y'] + ['x%d'%i for i in range(10)]]
    beta = randn(10)
    for t in range(1990, 2005):
        for r in c4:
            for c in c4[r]:
                x = [1] + list(randn(9))
                y = float(dot(beta, x) + randn(1))
                data.append([r, c, t, y] + list(x))

    fout = open('data.csv', 'w')
    csv.writer(fout).writerows(data)
    fout.close()
