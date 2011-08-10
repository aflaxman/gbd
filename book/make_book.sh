#!/bin/sh

python forward_sim.py
python forward_sim2.py
python rate_model.py
python poisson_model.py
python schiz_forest.py
python schiz_forest_plot.py

dexy
cd output
pdflatex theory-rate_model.tex
pdflatex theory-rate_model.tex

pdflatex theory-forward_sim.tex
pdflatex theory-forward_sim.tex

pdflatex ci-prev_meta_analysis.tex
pdflatex ci-prev_meta_analysis.tex
cd ..
