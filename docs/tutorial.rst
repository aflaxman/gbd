=====================
Quick-start Tutorials
=====================

Best-case scenario:  Lots of consistient, low-noise data
--------------------------------------------------------

`This Data Table <diabetes_data.tsv>`_ but it provides a
starting point for learning how to use DisMod.

1. Select all the rows and columns from the `simulated Type II Diabetes data table <diabetes_data.tsv>`_, and copy it (Edit->Copy or cntl-c).
2. Paste this into `Tab Separated Values` field on the `Data Upload Page </dismod/data/upload>`_.
3. Click `Create` to load this data into DisMod.  This will take you
   to the model summary page.  Try clicking on a cell in the sparkplot
   to the left of the table, to see a detailed view of some of the
   data you just loaded.
4. Click on `Adjust Priors` to set priors on the data.  In the `Adjust
   Priors` control panel:

   a. Select smoothness `moderately` for prevalence, incidence, and case fatality.
   b. Set zero range age before to 99 for remission, and set zero
      range age before to 1 for prevalence, incidence, and case fatality.
   c. Click `Apply`.

5. Click `Calculate/update covariates and age weights for model data`
   to load the age weights for the data set.  Note, this may take a
   while.

6. Click `Estimate empirical priors` to fit global pooled data for
   each parameter type separately.

7. Reload the ihme_dismod twitter page periodically.  When the IHME
   cluster finishes estimation for a portion of the model, it will
   produce a link directly to a summary graph.

8. Click `Run` on the sidebar to get back to the run-model
   page. Alternatively, use the browser back button.

9. Click `Fit each region-year-age individually` to
   generate a consistient estimate of disease parameters for each of the
   84 region/year/age triples in the GBD2005 study.

10. As before, reload the ihme_dismod twitter page periodically.  When the IHME
   cluster finishes estimation for a portion of the model, it will
   produce a link directly to a summary graph.
