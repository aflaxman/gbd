#!/bin/sh

python forward_sim.py
python forward_sim2.py

python binomial_model.py
python beta_binomial_model.py
python poisson_model.py
python neg_binom_sim.py
python schiz_forest.py

python age_patterns.py
python age_pattern_covariance.py
python age_pattern_mesh.py
python age_std.py

python af_age_groups.py
python bipolar.py

python hep_c.py
python hep_c_consistent.py
python hep_c_heterogeneity.py
python hep_c_smoothing.py

python pms.py
python pms_grics.py
python pms_ppc.py

dexy
cd output
pdflatex book.tex
pdflatex book.tex
cd ..
