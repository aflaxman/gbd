""" Covariate models"""

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

sex_value = {'male': .5, 'total':0., 'female': -.5}

def mean_covariate_model(name, mu, input_data, parameters, model, root_area, root_sex, root_year):
    """ Generate PyMC objects covariate adjusted version of mu

    Parameters
    ----------
    name : str
    mu : the unadjusted mean parameter for this node
    model : ModelData to use for covariates
    root_area, root_sex, root_year : str, str, int

    Results
    -------
    Returns dict of PyMC objects, including 'pi', the covariate
    adjusted predicted values for the mu and X provided
    """
    n = len(input_data.index)

    # make U and alpha
    p_U = model.hierarchy.number_of_nodes()  # random effects for area
    U = pandas.DataFrame(pl.zeros((n, p_U)), columns=model.hierarchy.nodes(), index=input_data.index)
    for i, row in input_data.T.iteritems():
        if row['area'] not in model.hierarchy:
            print 'WARNING: "%s" not in model hierarchy, skipping random effects for this observation' % row['area']
            continue
        
        for level, node in enumerate(nx.shortest_path(model.hierarchy, 'all', input_data.ix[i, 'area'])):
            model.hierarchy.node[node]['level'] = level
            U.ix[i, node] = 1.

    U = U.select(lambda col: U[col].std() > 1.e-5, axis=1)  # drop constant columns

    U_shift = pandas.Series(0., index=U.columns)
    for level, node in enumerate(nx.shortest_path(model.hierarchy, 'all', root_area)):
        if node in U_shift:
            U_shift[node] = 1.
    U = U - U_shift

    sigma_alpha = []
    for i in range(5):  # max depth of hierarchy is 5
        effect = 'sigma_alpha_%s_%d'%(name,i)
        if 'random_effects' in parameters and effect in parameters['random_effects']:
            prior = parameters['random_effects'][effect]
            print 'using stored RE for', effect, prior 
            sigma_alpha.append(mc.TruncatedNormal(effect, prior['mu'], pl.maximum(prior['sigma'], .001)**-2,
                                                  min(prior['mu'], prior['lower']),
                                                  max(prior['mu'], prior['upper']),
                                                  value=prior['mu']))
        else:
            sigma_alpha.append(mc.TruncatedNormal(effect, .1, .125**-2, .05, 5., value=.1))

    def MyTruncatedNormal(name, mu, tau, a, b, value):
        """ Need to make my own, because PyMC has underflow error when
        truncation is not doing anything"""
        @mc.stochastic(name=name)
        def my_trunc_norm(value=value, mu=mu, tau=tau, a=a, b=b):
            if tau**-.5 < (b-a)/10.:
                return mc.normal_like(value, mu, tau)
            else:
                return mc.truncated_normal_like(value, mu, tau, a, b)
        return my_trunc_norm
    
    alpha = pl.array([])
    alpha_potentials = []
    if len(U.columns) > 0:
        tau_alpha_index = []
        for alpha_name in U.columns:
            tau_alpha_index.append(model.hierarchy.node[alpha_name]['level'])
        tau_alpha_index=pl.array(tau_alpha_index, dtype=int)

        tau_alpha_for_alpha = [sigma_alpha[i]**-2 for i in tau_alpha_index]

        alpha = []
        for i, tau_alpha_i in enumerate(tau_alpha_for_alpha):
            effect = 'alpha_%s_%s'%(name, U.columns[i])
            if 'random_effects' in parameters and U.columns[i] in parameters['random_effects']:
                prior = parameters['random_effects'][U.columns[i]]
                print 'using stored RE for', effect, prior
                if prior['dist'] == 'TruncatedNormal':
                    alpha.append(MyTruncatedNormal(effect, prior['mu'], pl.maximum(prior['sigma'], .001)**-2,
                                                   prior['lower'], prior['upper'], value=0.))
                elif prior['dist'] == 'Constant':
                    alpha.append(float(prior['mu']))
                else:
                    assert 'ERROR: prior distribution "%s" is not implemented' % prior['dist']
            else:
                alpha.append(MyTruncatedNormal(effect, 0, tau=tau_alpha_i, a=-5., b=5., value=0))

        # change one stoch from each set of siblings in area hierarchy to a 'sum to zero' deterministic
        for parent in model.hierarchy:
            node_names = model.hierarchy.successors(parent)
            nodes = [U.columns.indexMap[n] for n in node_names if n in U]
            if len(nodes) > 0:
                i = nodes[0]
                old_alpha_i = alpha[i]

                # do not change if prior for this node has dist='constant'
                if parameters.get('random_effects', {}).get(U.columns[i], {}).get('dist') == 'Constant':
                    continue

                alpha[i] = mc.Lambda('alpha_det_%s_%d'%(name, i),
                                            lambda other_alphas_at_this_level=[alpha[n] for n in nodes[1:]]: -pl.sum(other_alphas_at_this_level))

                if isinstance(old_alpha_i, mc.Stochastic):
                    @mc.potential(name='alpha_pot_%s_%s'%(name, U.columns[i]))
                    def alpha_potential(alpha=alpha[i], mu=old_alpha_i.parents['mu'], tau=old_alpha_i.parents['tau'],
                                        a=old_alpha_i.parents['a'], b=old_alpha_i.parents['b']):
                        return mc.truncated_normal_like(alpha, mu, tau, a, b)
                    alpha_potentials.append(alpha_potential)

    # make X and beta
    X = input_data.select(lambda col: col.startswith('x_'), axis=1)

    # add sex as a fixed effect (TODO: decide if this should be in data.py, when loading gbd model)
    X['x_sex'] = [sex_value[row['sex']] for i, row in input_data.T.iteritems()]

    beta = pl.array([])
    X_shift = pandas.Series(0., index=X.columns)
    if len(X.columns) > 0:
        # shift columns to have zero for root covariate
        output_template = model.output_template.groupby(['area', 'sex', 'year']).mean()
        covs = output_template.filter(list(X.columns) + ['pop'])
        if len(covs.columns) > 1:
            leaves = [n for n in nx.traversal.bfs_tree(model.hierarchy, root_area) if model.hierarchy.successors(n) == []]
            if len(leaves) == 0:
                # networkx returns an empty list when the bfs tree is a single node
                leaves = [root_area]

            if root_sex == 'total' and root_year == 'all':  # special case for all years and sexes
                covs = covs.delevel().drop(['year', 'sex'], axis=1).groupby('area').mean()
                leaf_covs = covs.ix[leaves]
            else:
                leaf_covs = covs.ix[[(l, root_sex, root_year) for l in leaves]]

            for cov in covs:
                if cov != 'pop':
                    X_shift[cov] = (leaf_covs[cov] * leaf_covs['pop']).sum() / leaf_covs['pop'].sum()

        if 'x_sex' in X.columns:
            X_shift['x_sex'] = sex_value[root_sex]

        X = X - X_shift

        beta = []
        for i, effect in enumerate(X.columns):
            name_i = 'beta_%s_%d'%(name, i)
            if 'fixed_effects' in parameters and effect in parameters['fixed_effects']:
                prior = parameters['fixed_effects'][effect]
                print 'using stored FE for', effect, prior
                if prior['dist'] == 'normal':
                    beta.append(mc.Normal(name_i, mu=float(prior['mu']), tau=pl.maximum(prior['sigma'], .001)**-2, value=float(prior['mu'])))
                elif prior['dist'] == 'Constant':
                    beta.append(float(prior['mu']))
                else:
                    assert 'ERROR: prior distribution "%s" is not implemented' % prior['dist']
            else:
                beta.append(mc.Normal(name_i, mu=0., tau=.125**-2, value=0))

    @mc.deterministic(name='pi_%s'%name)
    def pi(mu=mu, U=pl.array(U, dtype=float), alpha=alpha, X=pl.array(X, dtype=float), beta=beta):
        return mu * pl.exp(pl.dot(U, alpha) + pl.dot(X, beta))

    return dict(pi=pi, U=U, U_shift=U_shift, sigma_alpha=sigma_alpha, alpha=alpha, alpha_potentials=alpha_potentials, X=X, X_shift=X_shift, beta=beta)



def dispersion_covariate_model(name, input_data, delta_lb, delta_ub):
    lower = pl.log(delta_lb)
    upper = pl.log(delta_ub)
    eta=mc.Uniform('eta_%s'%name, lower=lower, upper=upper, value=.5*(lower+upper))

    Z = input_data.select(lambda col: col.startswith('z_'), axis=1)
    Z = Z.select(lambda col: Z[col].std() > 0, 1)  # drop blank cols
    if len(Z.columns) > 0:
        zeta = mc.Normal('zeta', 0, .25**-2, value=pl.zeros(len(Z.columns)))

        @mc.deterministic(name='delta_%s'%name)
        def delta(eta=eta, zeta=zeta, Z=Z.__array__()):
            return pl.exp(eta + pl.dot(Z, zeta))

        return dict(eta=eta, Z=Z, zeta=zeta, delta=delta)

    else:
        @mc.deterministic(name='delta_%s'%name)
        def delta(eta=eta):
            return pl.exp(eta)
        return dict(eta=eta, delta=delta)



def predict_for(output_template, area_hierarchy, root_area, root_sex, root_year, area, sex, year, frac_unexplained, vars, lower, upper):
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
    frac_unexplained : float, fraction of unexplained variation to include in prediction
    vars : dict, including entries for alpha, beta, mu_age, U, and X
    lower, upper : float, bounds on predictions from expert priors

    Results
    -------
    Returns array of draws from posterior predicted distribution
    """
    if 'alpha' in vars and isinstance(vars['alpha'], mc.Node):
        alpha_trace = vars['alpha'].trace()
    elif 'alpha' in vars and isinstance(vars['alpha'], list):
        
        alpha_trace = []
        for n in vars['alpha']:
            if isinstance(n, mc.Stochastic):
                alpha_trace.append(n.trace())
            else:
                alpha_trace.append([float(n) for i in vars['mu_age'].trace()])
        alpha_trace = pl.vstack(alpha_trace).T
    else:
        alpha_trace = pl.array([])

    if 'beta' in vars and isinstance(vars['beta'], mc.Node):
        beta_trace = vars['beta'].trace()
    elif 'beta' in vars and isinstance(vars['beta'], list):
        
        beta_trace = []
        for n in vars['beta']:
            if isinstance(n, mc.Stochastic):
                beta_trace.append(n.trace())
            else:
                beta_trace.append([float(n) for i in vars['mu_age'].trace()])
        beta_trace = pl.vstack(beta_trace).T
    else:
        beta_trace = pl.array([])

    # commenting this out leads to including country random effect from sigma_alpha if there is not data
    #if len(alpha_trace) == 0 and len(beta_trace) == 0:
    #    # TODO: create country level random effects here, based on sigma_alpha
    #    return vars['mu_age'].trace()

    leaves = [n for n in nx.traversal.bfs_tree(area_hierarchy, area) if area_hierarchy.successors(n) == []]
    if len(leaves) == 0:
        # networkx returns an empty list when the bfs tree is a single node
        leaves = [area]

    covariate_shift = 0.
    total_population = 0.

    output_template = output_template.groupby(['area', 'sex', 'year']).mean()
    covs = output_template.filter(vars['X'].columns)
    if 'x_sex' in vars['X'].columns:
        covs['x_sex'] = sex_value[sex]
    assert pl.all(covs.columns == vars['X_shift'].index), 'covariate columns and unshift index should match up'
    for x_i in vars['X_shift'].index:
        covs[x_i] -= vars['X_shift'][x_i] # shift covariates so that the root node has X_ar,sr,yr == 0
    

    # make U_l, outside of loop, but initialize inside loop
    # this allows same random effect draws across countries; necessary?
    p_U = area_hierarchy.number_of_nodes()  # random effects for area
    U_l = pandas.DataFrame(pl.zeros((1, p_U)), columns=area_hierarchy.nodes())
    U_l = U_l.filter(vars['U'].columns)

    for l in leaves:
        log_shift_l = 0.
        U_l.ix[0,:] = 0.
        
        for node in nx.shortest_path(area_hierarchy, root_area, l):
            if node not in U_l.columns:
                ## Add a columns U_l[node] = rnormal(0, appropriate_tau)
                level = len(nx.shortest_path(area_hierarchy, 'all', node))-1
                tau_l = vars['sigma_alpha'][level].trace()**-2
                U_l[node] = 0.
                if len(alpha_trace) > 0:
                    alpha_trace = pl.vstack((alpha_trace.T, mc.rnormal(0., tau_l))).T
                else:
                    alpha_trace = mc.rnormal(0., tau_l)
            U_l.ix[0, node] = 1.

        for node in vars['U_shift']:
            U_l -= vars['U_shift'][node]
        
        log_shift_l += pl.dot(alpha_trace, pl.atleast_2d(U_l).T)
            
        # make X_l
        if len(beta_trace) > 0:
            X_l = covs.ix[l, sex, year]
            log_shift_l += pl.dot(beta_trace, pl.atleast_2d(X_l).T)

        shift_l = pl.exp(log_shift_l)
        covariate_shift += shift_l * output_template['pop'][l,sex,year]
        total_population += output_template['pop'][l,sex,year]
    covariate_shift /= total_population

    parameter_prediction = vars['mu_age'].trace() * covariate_shift

    # add requested portion of unexplained variation
    additional_shift = 0.

    ## uncomment below to include additional variation from negative binomial overdispersion
    #if 'eta' in vars and frac_unexplained > 0.:
    #    delta_trace = pl.exp(vars['eta'].trace())
    #    var = pl.log((1+delta_trace) / delta_trace)
    #    additional_shift = mc.rnormal(0., 1. / (var * frac_unexplained))

    parameter_prediction = (parameter_prediction.T * pl.exp(additional_shift)).T

    # clip predictions to bounds from expert priors
    parameter_prediction = parameter_prediction.clip(lower, upper)
    
    return parameter_prediction
    
