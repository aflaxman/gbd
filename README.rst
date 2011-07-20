Introduction
============

This special, stripped-down version of DisMod III is designed to fit
models with birth prevalence, which is known (mean and uncertainty).

In this simplified model, there is no need to borrow strength between
regions, years, or sexes, and little more than forward simulation is
needed, starting with birth prevalence levels dictated by a single
input data point.

This means the model can run much faster, and provides as opportunity
to strip out the unneeded code and carefully check the code that
remains.

