

# with no excess-mortality or remission
dm = dismod3.load_disease_model(20088)
fit_posterior.fit_posterior(dm, 'asia_east', 'female', '2005', store_results=False)
book_graphics.plot_age_patterns(dm, region='asia_east', sex='female', rate_types=['incidence', 'prevalence'])

# with default excess-mortality and remission
dm = dismod3.load_disease_model(20087)
fit_posterior.fit_posterior(dm, 'asia_east', 'female', '2005', store_results=False)
book_graphics.plot_age_patterns(dm, region='asia_east', sex='female', rate_types=['incidence', 'prevalence'])
