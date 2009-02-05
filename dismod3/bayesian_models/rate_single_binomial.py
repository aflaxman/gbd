from dismod3.models.probabilistic_utils import *

def rate_model_vars(asrf):
    vars = {}

    add_rate_stochs(vars, 'asrf', asrf.fit['age_mesh'], asrf.fit['out_age_mesh'])

    @mc.potential
    def smoothing_logp(f=vars['logit(asrf)'], tau=1./(.5)**2):
        return mc.normal_like(np.diff(f), 0.0, tau)

    @mc.potential
    def zero_logp(f=vars['asrf'], age_start=0, age_end=5):
        return mc.normal_like(f[range(age_start, age_end)], 0.0, 1./(1e-3)**2)

    vars['priors'] = [smoothing_logp, zero_logp]

    vars['observed_rates'] = observed_rates_stochs(asrf.rates.all(), vars['asrf'])

    return vars

