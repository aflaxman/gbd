===================================================
Introduction to the Generic Disease Modeling System
===================================================

DisMod III is IHME's newest iteration of the generic disease modeling
system. It has been redesigned from the ground up using a fully
Bayesian model.  DisMod III exists to produce consistent estimates for
epidemiological parameters of disease.  It will be used to provide
estimates of age-specific incidence and duration of the 200+ diseases,
injuries, and risk factors that are being included in the Global
Burden of Disease 2005 (GBD2005) Study; in GBD2005, incidence,
duration, and disability weights will be the inputs for calculating
Years Lived with Disability in the Disability Adjusted Life Year
measurement of disease burden.  DisMod III has been designed to be as
friendly and accurate as possible.

--------------------
Why a Generic Model?
--------------------

GBD2005 studies over 200 different diseases, injuries, and risk
factors for 22 regions of the world.  Even for well-known diseases in
developed regions available epidemiological data is sparse and noisy.
By applying a simple four-compartment model of how disease moves
through a population, researchers can combine data from multiple
sources to reduce errors in the data that is available, reconcile the
data that is inconsistent, and impute or forecast the data that is
absent altogether.

------------------
DisMod is Friendly
------------------

DisMod III runs as a web-based application, allowing disease experts
to analyze descriptive epidemiological data in their web-browsers,
without installing any additional software.  It has been developed as
Free/Libre Open-Source Software (FLOSS), to permit external code
review and research replication.  DisMod III is designed to be used.
Initially it will be used internally by IHME researchers, and as it is
developed, it is intended to be useful more widely for
epidemiologists, doctors, and public health researchers.  It will be
available from IHME as a web-service, and also available as a
Python/Django application for researchers who prefer to host their own
installation.  DisMod III follows an Open Research pattern, making it
easy to record and annotate the process of exploratory data analysis,
model selection, and prior specification.

------------------
DisMod is Accurate
------------------

All models are wrong, but some are worse than useless, by providing
false confidence about wrong results.  DisMod III will strive to
provide estimates of epidemiological parameters with accurate
uncertainty intervals.  It is important to generate point estimates of
parameters which are valid, but it is also important to generate
uncertainty intervals that are not more certain than they should be.
DisMod III uses a beta-binomial model of disease prevalence,
incidence, remission, and case-fatality rates to attempt to avoid the
erroneously small uncertainty intervals that could be caused by model
misspecification.  Because DisMod III fits models with a randomized
Markov-Chain Monte Carlo algorithm (using the Adaptive Metropolis step
function), it is important to include convergence diagnostics and
goodness-of-fit tests, and to validate DisMod III output extensively
against realistic examples that have been simulated in a setting where
the estimates can be compared to ground truth.
