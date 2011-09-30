""" Covariate models"""

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

sex_value = {'male': 1., 'total':0., 'female': -1.}

def mean_covariate_model(name, mu, data, output_template, area_hierarchy, root_area, root_sex, root_year):
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

    # make U and alpha
    p_U = 2 + area_hierarchy.number_of_nodes()  # random effects for sex, time, area
    U = pandas.DataFrame(pl.zeros((n, p_U)), columns=['sex', 'time'] + area_hierarchy.nodes(), index=data.index)
    for i, row in data.T.iteritems():
        U.ix[i, 'sex'] = sex_value[data.ix[i, 'sex']] - sex_value[root_sex]
        U.ix[i, 'time'] = .5 * (data.ix[i, 'year_start'] + data.ix[i, 'year_end'])
        if root_year == 'all':
             U.ix[i, 'time'] -= 2000.
        else:
            U.ix[i, 'time'] -= root_year
        for node in nx.shortest_path(area_hierarchy, root_area, data.ix[i, 'area']):
            # do not include random effect for lowest level of hierarchy
            if len(area_hierarchy.successors(node)) == 0:
                continue

            U.ix[i, node] = 1.
    U = U.select(lambda col: U[col].std() > 1.e-5, axis=1)  # drop blank columns

    tau_alpha = pl.array([])
    alpha = pl.array([])
    if len(U.columns) > 0:
        tau_alpha = mc.InverseGamma(name='tau_alpha_%s'%name, alpha=1., beta=1., value=pl.ones_like(U.columns))
        # TODO: consider parameterizations where tau_alpha is the same for different areas (or just for different areas that are children of the same area)
        alpha = mc.Normal(name='alpha_%s'%name, mu=0, tau=.01**-2, value=pl.zeros(len(U.columns)))  

    # make X and beta
    X = data.select(lambda col: col.startswith('x_'), axis=1)
    X = X.select(lambda col: X[col].std() > 1.e-5, axis=1)  # drop blank columns

    beta = pl.array([])
    X_shift = 0
    if len(X.columns) > 0:
        # shift columns to have zero for root covariate
        output_template = output_template.groupby(['area', 'sex', 'year']).mean()
        covs = output_template.filter(X.columns)
        if len(covs.columns) > 0:
            leaves = [n for n in nx.traversal.bfs_tree(area_hierarchy, root_area) if area_hierarchy.successors(n) == []]
            if len(leaves) == 0:
                # networkx returns an empty list when the bfs tree is a single node
                leaves = [root_area]
            if root_sex == 'total' and root_year == 'all':  # special case for all years and sexes
                covs = covs.delevel().drop(['year', 'sex'], axis=1).groupby('area').mean()
                X_shift = covs.ix[leaves].mean()
            else:
                X_shift = covs.ix[[(l, root_sex, root_year) for l in leaves]].mean()

            X = X - X_shift

        #beta = mc.Uniform('beta_%s'%name, -5., 5., value=pl.zeros(len(X.columns)))
        beta = mc.Normal('beta_%s'%name, mu=0., tau=.01**-2, value=pl.zeros(len(X.columns)))

    @mc.deterministic(name='pi_%s'%name)
    def pi(mu=mu, U=pl.array(U, dtype=float), alpha=alpha, X=pl.array(X, dtype=float), beta=beta):
        return mu * pl.exp(pl.dot(U, alpha) + pl.dot(X, beta))

    return dict(pi=pi, U=U, tau_alpha=tau_alpha, alpha=alpha, X=X, X_shift=X_shift, beta=beta)

def dispersion_covariate_model(name, data):

    eta=mc.Normal('eta_%s'%name, mu=5., tau=1., value=5.)

    Z = data.select(lambda col: col.startswith('z_'), axis=1)
    Z = Z.select(lambda col: Z[col].std() > 0, 1)  # drop blank cols
    if len(Z.columns) > 0:
        zeta = mc.Uniform('zeta', -5, 5, value=pl.zeros(len(Z.columns)))

        @mc.deterministic(name='delta_%s'%name)
        def delta(eta=eta, zeta=zeta, Z=Z.__array__()):
            return (50. + pl.exp(eta)) * pl.exp(pl.dot(Z, zeta))

        return dict(eta=eta, Z=Z, zeta=zeta, delta=delta)

    else:
        @mc.deterministic(name='delta_%s'%name)
        def delta(eta=eta):
            return (50. + pl.exp(eta))
        return dict(eta=eta, delta=delta)



def predict_for(output_template, area_hierarchy, root_area, root_sex, root_year, area, sex, year, vars):
    """ Generate draws from posterior predicted distribution for a
    specific (area, sex, year)

    Parameters
    ----------
    output_template : pandas.DataFrame with covariate data for all leaf nodes in area hierarchy
    area_hierarchy : nx.DiGraph encoding hierarchical relationship of areas
    root_area : str, area for which this model was fit consistently
    root_sex : str, area for which this model was fit consistently
    root_year : str, area for which this model was fit consistently
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

    leaves = [n for n in nx.traversal.bfs_tree(area_hierarchy, area) if area_hierarchy.successors(n) == []]
    if len(leaves) == 0:
        # networkx returns an empty list when the bfs tree is a single node
        leaves = [area]

    covariate_shift = 0.
    total_population = 0.

    p_U = 2 + area_hierarchy.number_of_nodes()  # random effects for sex, time, area
    U_l = pandas.DataFrame(pl.zeros((1, p_U)), columns=['sex', 'time'] + area_hierarchy.nodes())
    U_l = U_l.filter(vars['U'].columns)
    
    output_template = output_template.groupby(['area', 'sex', 'year']).mean()
    covs = output_template.filter(vars['X'].columns)
    covs -= vars['X_shift'] # shift covariates so that the root node has X_ar,sr,yr == 0
    
    for l in leaves:
        log_shift_l = 0.

        # make U_l
        if len(alpha_trace) > 0:
            U_l.ix[0, :] = 0.
            if 'sex' in U_l:
                U_l.ix[0, 'sex'] = sex_value[sex] - sex_value[root_sex]
                
            if root_year == 'all':
                U_l.ix[0, 'time'] = year - 2000.
            else:
                U_l.ix[0, 'time'] = year - root_year

            for node in nx.shortest_path(area_hierarchy, root_area, l):
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
    
