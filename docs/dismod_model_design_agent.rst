=========================
DisMod Model Design Agent
=========================

DisMod has reached the point where it is mature enough to rethink how
analysts input their data, select the covariates included in the
model, and specify their expert priors.

The goal of this document is to precisely specify the requirements for
our new DM-MDA.  It is intended to be an collection of incremental
changes to the existing java application for data checking and java
applets for expert prior elicitation and covariate selection.


Reorganization of the Expert Priors
-----------------------------------

1.  The expert priors are now organized by prior type and then by
    parameter type; I think reversing them makes sense: choose a parameter
    type, for example remission, and then work on setting the priors for
    that parameter.

2.  Since merging country-level covariates is now handled in the data
    checker, the covariate selection applet can be simplified.  It does
    not need to distinguish between country-level covariates and
    study-level covariates, and it should be specific to each parameter
    type.  For example, it should be possible to set prevalence to use
    LDI_pc as a covariate without using it for remission.

3.  The data and parameters selected by the user must be saved in a
    standard format, which must be fully specified.  It should bear a
    similarity to the disease_json format already specified, with a list
    of data dicts and a parameter dict, but it can now be more precise.

4.  The DM-DMA should be able to load a disease model that is saved in
    the format from section 3, for additional parameter or data
    modifications.

5.  The DM parameters and priors should be compatible with the
    parameters needed by Brad's DM4 model, but appropriately simpler.

6.  The DM-DMA should allow reorganizing the hierarchical structure of
    the "atomic units" of the analysis, i.e. changing the structure of the
    regions and super-regions.  This hierarchy should be stored as part of
    the disease format specified in requirement 3

7.  The DM-DMA should prepare an "model output data template" in a
    format analogous to the dismod data csv, which specifies which
    country/sex/year/age ranges should be produced as output of the model,
    and the country-specific covariates necessary for generating these
    estimates should be merged into the model at the same time as the
    country-level covariates are merged into the input data.
