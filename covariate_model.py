""" Covariate models"""

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

sex_value = {'male': 1., 'total':0., 'female': -1.}

def mean_covariate_model(name, mu, data, hierarchy, root):
    """ Generate PyMC objects covariate adjusted version of mu

    Parameters
    ----------
    name : str
    mu : the unadjusted mean parameter for this node
    beta : effect coefficient prior
    data : pandas.DataFrame containing design matrix for covariates,
    as well as columns for year_start, year_end, and area
    hierarchy : nx.DiGraph encoding hierarchical structure and similarity weights

    Results
    -------
    Returns dict of PyMC objects, including 'pi', the covariate
    adjusted predicted values for the mu and X provided

    Notes
    -----
    This is used twice, first to predict the means of the observed
    data, and then to predict the level values for the node
    """
    n = len(data.index)

    # make U
    p_U = 2 + hierarchy.number_of_nodes()  # random effects for sex, time, area
    U = pandas.DataFrame(pl.zeros((n, p_U)), columns=['sex', 'time'] + hierarchy.nodes(), index=data.index)
    for i, row in data.T.iteritems():
        U.ix[i, 'sex'] = sex_value[data.ix[i, 'sex']]
        U.ix[i, 'time'] = .5 * (data.ix[i, 'year_start'] + data.ix[i, 'year_end']) - 2000.
        for node in nx.shortest_path(hierarchy, root, data.ix[i, 'area']):
            U.ix[i, node] = 1.
    U = U.select(lambda col: U[col].std() > 1.e-5, axis=1)  # drop blank columns

    # shift columns to have mean zero?
    #U_shift = U.mean(0)
    #U = U - U_shift

    # make tau_alpha and alpha
    if len(U.columns) > 0:
        tau_alpha = mc.InverseGamma(name='tau_alpha_%s'%name, alpha=1., beta=1., value=pl.ones_like(U.columns))  # TODO: consider parameterizations where tau_alpha is the same for different areas (or just for different areas that are children of the same area)
        alpha = mc.Normal(name='alpha_%s'%name, mu=0, tau=.01**-2, value=pl.zeros(len(U.columns)))  
    else:
        tau_alpha = pl.array([])
        alpha = pl.array([])

    # TODO: consider faster ways to calculate dot(U, alpha), since the matrix is sparse and (half-)integral

    X = data.select(lambda col: col.startswith('x_'), axis=1)
    X = X.select(lambda col: X[col].std() > 1.e-5, axis=1)  # drop blank columns
    # shift columns to have mean zero?
    X_shift = X.mean(0)
    X = X - X_shift
    if len(X.columns) > 0:
        #beta = mc.Uniform('beta_%s'%name, -5., 5., value=pl.zeros(len(X.columns)))
        beta = mc.Normal('beta_%s'%name, mu=0., tau=.05**-2, value=pl.zeros(len(X.columns)))
    else:
        beta = pl.array([])

    @mc.deterministic(name='pi_%s'%name)
    def pi(mu=mu, U=pl.array(U, dtype=float), alpha=alpha, X=pl.array(X, dtype=float), beta=beta):
        return mu * pl.exp(pl.dot(U, alpha) + pl.dot(X, beta))

    return dict(pi=pi, U=U, tau_alpha=tau_alpha, alpha=alpha, X=X, beta=beta)

def dispersion_covariate_model(name, data):

    eta=mc.Normal('eta_%s'%name, mu=5., tau=1., value=5.)

    Z = data.select(lambda col: col.startswith('z_'), axis=1)
    Z = Z.select(lambda col: Z[col].std() > 0, 1)  # drop blank cols
    if len(Z.columns) > 0:
        zeta = mc.Uniform('zeta', -5, 5, value=pl.zeros(len(Z.columns)))
    else:
        zeta = pl.array([])

    @mc.deterministic(name='delta_%s'%name)
    def delta(eta=eta, zeta=zeta, Z=Z.__array__()):
        return (50. + pl.exp(eta)) * pl.exp(pl.dot(Z, zeta))

    return dict(eta=eta, Z=Z, zeta=zeta, delta=delta)


def predict_for(output_template, hierarchy, root, area, sex, year, vars):
    """ Generate draws from posterior predicted distribution for a
    specific (area, sex, year)

    Parameters
    ----------
    output_template : pandas.DataFrame with covariate data for all leaf nodes in area hierarchy
    hierarchy : nx.DiGraph encoding hierarchical relationship of areas
    root : str, area for which this model was fit consistently
    area : str, area to predict for
    sex : str, sex to predict for
    year : str, year to predict for
    vars : dict, including entries for alpha, beta, mu_age, U, and X

    Results
    -------
    Returns array of draws from posterior predicted distribution
    """

    if 'alpha' in vars and isinstance(vars['alpha'], mc.Node):
        alpha_trace = vars['alpha'].trace()
    else:
        alpha_trace = pl.array([])

    if 'beta' in vars and isinstance(vars['beta'], mc.Node):
        beta_trace = vars['beta'].trace()
    else:
        beta_trace = pl.array([])

    if len(alpha_trace) == 0 and len(beta_trace) == 0:
        return vars['mu_age'].trace()

    leaves = [n for n in nx.traversal.bfs_tree(hierarchy, area) if hierarchy.successors(n) == []]
    if len(leaves) == 0:
        # networkx returns an empty list when the bfs tree is a single node
        leaves = [area]

    covariate_shift = 0.
    total_population = 0.

    p_U = 2 + hierarchy.number_of_nodes()  # random effects for sex, time, area
    U_l = pandas.DataFrame(pl.zeros((1, p_U)), columns=['sex', 'time'] + hierarchy.nodes())
    U_l = U_l.filter(vars['U'].columns)
    
    output_template = output_template.groupby(['area', 'sex', 'year']).mean()
    covs = output_template.filter(vars['X'].columns)
    
    for l in leaves:
        log_shift_l = 0.

        # make U_l
        if len(alpha_trace) > 0:
            U_l.ix[0, :] = 0.
            if 'sex' in U_l:
                U_l.ix[0, 'sex'] = covariate_model.sex_value[sex]

            U_l.ix[0, 'time'] = year - 2000.
            for node in nx.shortest_path(hierarchy, root, l):
                if node in U_l.columns:
                    U_l.ix[0, node] = 1.
                else:
                    # TODO: include appropriate uncertainty for random effects that are not in model
                    pass

            log_shift_l += pl.dot(alpha_trace, pl.atleast_2d(U_l).T)
            
        # make X_l
        if len(beta_trace) > 0:
            X_l = covs.ix[l, sex, year]
            log_shift_l += pl.dot(beta_trace, pl.atleast_2d(X_l).T)

        shift_l = pl.exp(log_shift_l)
        covariate_shift += shift_l * output_template['pop'][l,sex,year]
        total_population += output_template['pop'][l,sex,year]
    covariate_shift /= total_population

    return vars['mu_age'].trace() * covariate_shift
    
