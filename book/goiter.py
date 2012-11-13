import dismod3

dismod3.data.fetch_disease_model_if_necessary(29824, 'models/goiter')

model = dismod3.data.load('models/goiter/')
model.keep(areas=['IRN'])
model.vars += dismod3.ism.age_specific_rate(model, 'p')
dismod3.fit.fit_asr(model)
graphics.plot_one_type(model, model.vars['p'], {}, 'p')
model.vars.plot_acorr()


model = dismod3.data.load('models/goiter/')
model.keep(areas=['north_africa_middle_east'])
model.vars += dismod3.ism.age_specific_rate(model, 'p')
dismod3.fit.fit_asr(model)
graphics.plot_one_type(model, model.vars['p'], {}, 'p')
model.vars.plot_acorr()


model = dismod3.data.load('models/goiter/')
model.keep(areas=['super-region_%d'%i for i in range(2,8)])
model.vars += dismod3.ism.age_specific_rate(model, 'p')
dismod3.fit.fit_asr(model)
graphics.plot_one_type(model, model.vars['p'], {}, 'p')
model.vars.plot_acorr()


model = dismod3.data.load('models/goiter/')
model.keep(areas=['super-region_1'])
model.vars += dismod3.ism.age_specific_rate(model, 'p')
dismod3.fit.fit_asr(model)
model.vars.plot_acorr()
graphics.plot_one_type(model, model.vars['p'], {}, 'p')


model = dismod3.data.load('models/goiter/')
model.keep(areas=['super-region_0'])
model.vars += dismod3.ism.age_specific_rate(model, 'p')
dismod3.fit.fit_asr(model)
model.vars.plot_acorr()
graphics.plot_one_type(model, model.vars['p'], {}, 'p')


model = dismod3.data.load('models/goiter/')
model.vars += dismod3.ism.age_specific_rate(model, 'p')
dismod3.fit.fit_asr(model)
model.vars.plot_acorr()
graphics.plot_one_type(model, model.vars['p'], {}, 'p')
raphics.plot_one_effects(model.vars['p'], 'p', model.hierarchy)
