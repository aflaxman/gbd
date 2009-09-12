=====================
Quick-start Tutorials
=====================

Best-case scenario:  Lots of consistient, low-noise data
--------------------------------------------------------

`This Data Table <diabetes_data.tsv>`_ provides a starting point for
learning how to use DisMod.  It is a simple example of synthetically
generated data inspired by actual Diabetes incidence and relative risk
data from Thailand.

1. Select all the rows and columns from the `simulated Type II Diabetes data table <diabetes_data.tsv>`_, and copy it (Edit->Copy or cntl-c).
2. Paste this into `Tab Separated Values` field on the `Data Upload Page </dismod/data/upload>`_.
3. Click `Create` to load this data into DisMod.  This will take you
   to the model summary page.  The specks of color to the left of the
   table are a tiny summary of all the data you just loaded, in a format that
   Tufte might call the small-multiples sparkplot.  Try clicking on a cell in
   the sparkplot to see a detailed view of some of the
   data you just loaded.
4. Click on `Adjust Priors` to set priors on the data.  In the `Adjust
   Priors` control panel:

   a. Select smoothness `moderately` for prevalence, incidence, and case fatality.
   b. Set zero range age before to 99 for remission, and set zero
      range age before to 1 for prevalence, incidence, and case fatality.
   c. Click `Apply`.

5. Click `Calculate covariates for model data`
   to load the age weights for the data set.  Note, this may take a
   while.

6. Click `Estimate empirical priors` to fit global pooled data for
   each parameter type separately.

7. Be patient; it should take about 1 minute for DisMod to fit the
   empirical priors for this model.  As DisMod works, this webpage
   will refresh every 30 seconds to show the results so far.  Note:
   DisMod fits empirical priors for incidence, prevalence, and case
   fatality in parallel, and does not yet have a mechanism to learn when
   all these computations are complete.  That is why is says "queued",
   "running", and "(at least partially) finished"...  DisMod doesn't know
   when it is all done, but _you_ do.

9. Click `Estimate posterior` to
   generate a consistient estimate of disease parameters for some of the
   84 region/year/age triples in the GBD2005 study.  (Note: in the
   future you will be able to select specific regions to fit, but for now
   it will fit only Southeast Asia, which is the only region with data
   in this tutorial.  Also to speed things up, DisMod will not estimate
   the uncertainty interval around this posterior point estimate right
   now.  But rest assured that it is capable of doing so.).

10. As before, be patient.  For this detailed computation, DisMod can take up to 20 minutes.
