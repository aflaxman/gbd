#!/bin/sh

time python forward_sim.py
time python forward_sim2.py

time python binomial_model.py
time python beta_binomial_model.py
time python poisson_model.py
time python neg_binom_sim.py
time python schiz_forest.py
time python zero_forest.py

time python age_patterns.py
time python age_pattern_covariance.py
time python age_pattern_mesh.py
time python age_std.py

time python af_age_groups.py
time python bipolar.py

time python hep_c.py
time python hep_c_consistent.py
time python hep_c_heterogeneity.py
time python hep_c_smoothing.py

time python pms.py
time python pms_grids.py
time python pms_ppc.py

time dexy
cd output
pdflatex book.tex
pdflatex book.tex
cd ..
