""" Stub module for planned refactoring of dismod3 model creation methods"""


# TODO: refactor data_model.data_model into ism.age_specific_rate
import data_model
import data
def age_specific_rate(model, data_type, reference_area='all', reference_sex='total', reference_year='all'):
    result = data.ModelVars()
    result[data_type] = data_model.data_model(data_type, model, data_type,
                                              reference_area, reference_sex, reference_year,
                                              None, None, None)
    return result


# TODO: refactor consistent_model.consistent_model into ism.consistent
import consistent_model
reload(consistent_model)
def consistent(model, reference_area='all', reference_sex='total', reference_year='all'):
    priors = {}
    return consistent_model.consistent_model(model, reference_area, reference_sex, reference_year, priors)
