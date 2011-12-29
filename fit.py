""" Stub module for planned refactoring of dismod3 model fitting methods"""

# TODO: move fit_model.fit_data_model to fit.fit_asr
import fit_model
def fit_asr(model, iter=20000, burn=10000, thin=10, tune_interval=1000):
    for k in model.vars:
        if k != 'logit_C0':
            fit_model.fit_data_model(model.vars[k], iter, burn, thin, tune_interval)
