""" Stub module for planned refactoring of dismod3 model fitting methods"""

# TODO: move fit_model.fit_data_model to fit.fit_asr
import fit_model
def fit_asr(vars, iter, burn, thin, tune_interval):
    for k in vars:
        if k != 'logit_C0':
            fit_model.fit_data_model(vars[k], iter, burn, thin, tune_interval)
