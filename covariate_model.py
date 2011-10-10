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
    p_U = 1 + area_hierarchy.number_of_nodes()  # random effects for sex, time, area
    U = pandas.DataFrame(pl.zeros((n, p_U)), columns=['sex'] + area_hierarchy.nodes(), index=data.index)
    for i, row in data.T.iteritems():
        U.ix[i, 'sex'] = sex_value[data.ix[i, 'sex']] - sex_value[root_sex]
        for level, node in enumerate(nx.shortest_path(area_hierarchy, root_area, data.ix[i, 'area'])):
            area_hierarchy.node[node]['level'] = level
            U.ix[i, node] = 1.

    U = U.select(lambda col: U[col].std() > 1.e-5, axis=1)  # drop constant columns

    sigma_alpha = [mc.TruncatedNormal('sigma_alpha_%s_%d'%(name,i), .003, .125**-2, .001, .25, value=.003) for i in range(5)]  # max depth of hierarchy is 4
    alpha = pl.array([])
    alpha_potentials = []
    if len(U.columns) > 0:
        tau_alpha_index = []
        for alpha_name in U.columns:
            if alpha_name == 'sex':
                tau_alpha_index.append(0)
            else:
                tau_alpha_index.append(area_hierarchy.node[alpha_name]['level']+1)
        tau_alpha_index=pl.array(tau_alpha_index, dtype=int)

        tau_alpha_for_alpha = [sigma_alpha[i]**-2 for i in tau_alpha_index]
        alpha = [mc.TruncatedNormal(name='alpha_%s_%d'%(name, i), mu=0, tau=tau_alpha_i, a=-.5, b=.5, value=0) for i, tau_alpha_i in enumerate(tau_alpha_for_alpha)]

        # change one stoch from each set of siblings in area hierarchy to a 'sum to zero' deterministic
        for parent in area_hierarchy:
            node_names = area_hierarchy.successors(parent)
            nodes = [U.columns.indexMap[n] for n in node_names if n in U]
            if len(nodes) > 0:
                i = nodes[0]
                alpha[i] = mc.Lambda('alpha_det_%s_%d'%(name, i),
                                            lambda other_alphas_at_this_level=[alpha[n] for n in nodes[1:]]: -pl.sum(other_alphas_at_this_level))

                @mc.potential(name='alpha_pot_%s_%d'%(name, i))
                def alpha_potential(alpha=alpha[nodes[0]], tau=tau_alpha_for_alpha[i]):
                    return mc.truncated_normal_like(alpha, 0, tau, -.5, .5) 
                alpha_potentials.append(alpha_potential)

    # make X and beta
    X = data.select(lambda col: col.startswith('x_'), axis=1)
    X = X.select(lambda col: X[col].std() > 1.e-5, axis=1)  # drop blank columns

    beta = pl.array([])
    X_shift = pandas.DataFrame()
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

        beta = [mc.Normal('beta_%s_%d'%(name, i), mu=0., tau=.125**-2, value=0) for i in range(len(X.columns))]

    @mc.deterministic(name='pi_%s'%name)
    def pi(mu=mu, U=pl.array(U, dtype=float), alpha=alpha, X=pl.array(X, dtype=float), beta=beta):
        return mu * pl.exp(pl.dot(U, alpha) + pl.dot(X, beta))

    return dict(pi=pi, U=U, sigma_alpha=sigma_alpha, alpha=alpha, alpha_potentials=alpha_potentials, X=X, X_shift=X_shift, beta=beta)



def dispersion_covariate_model(name, data):

    eta=mc.Normal('eta_%s'%name, mu=5., tau=.25**-2, value=5.)

    Z = data.select(lambda col: col.startswith('z_'), axis=1)
    Z = Z.select(lambda col: Z[col].std() > 0, 1)  # drop blank cols
    if len(Z.columns) > 0:
        zeta = mc.Normal('zeta', 0, .25**-2, value=pl.zeros(len(Z.columns)))

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
    elif 'alpha' in vars and isinstance(vars['alpha'], list):
        alpha_trace = pl.vstack([n.trace() for n in vars['alpha']]).T
    else:
        alpha_trace = pl.array([])

    if 'beta' in vars and isinstance(vars['beta'], mc.Node):
        beta_trace = vars['beta'].trace()
    elif 'beta' in vars and isinstance(vars['beta'], list):
        beta_trace = pl.vstack([n.trace() for n in vars['beta']]).T
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
                
            if 'time' in U_l:
                if root_year == 'all':
                    U_l.ix[0, 'time'] = year - 2000.
                else:
                    U_l.ix[0, 'time'] = year - root_year

            for node in nx.shortest_path(area_hierarchy, root_area, l):
                if node not in U_l.columns:
                    ## Add a columns U_l[node] = rnormal(0, appropriate_tau)
                    level = 1 + (len(nx.shortest_path(area_hierarchy, root_area, node))-1)
                    tau_l = vars['sigma_alpha'][level].trace()**-2
                    U_l[node] = 0.
                    alpha_trace = pl.vstack((alpha_trace.T, mc.rnormal(0., tau_l))).T
                U_l.ix[0, node] = 1.

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
    
