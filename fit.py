""" Stub module for planned refactoring of dismod3 model fitting methods"""

# TODO: move fit_model.fit_data_model to fit.fit_asr
import fit_model
reload(fit_model)
def fit_asr(model, data_type, iter=20000, burn=10000, thin=10, tune_interval=1000, verbose=True):
    model.map, model.mcmc = fit_model.fit_data_model(model.vars[data_type], iter, burn, thin, tune_interval, verbose)

# TODO: move fit_model.fit_consistent_model to fit.fit_consistent
import fit_model
def fit_consistent(model, iter=20000, burn=10000, thin=10, tune_interval=1000, verbose=True):
    model.map, model.mcmc = fit_model.fit_consistent_model(model.vars, iter, burn, thin, tune_interval, verbose)

