=====================
Quick-start Tutorials
=====================

Best-case scenario:  Lots of consistient, low-noise data
--------------------------------------------------------

`This Data Table <dismoditis.tsv>`_ provides a starting point for
learning how to use DisMod.  It is a simple example of synthetically
generated data inspired by Type II Diabetes incidence and relative risk
data from Thailand.

1. Select all the rows and columns from the `data table <dismoditis.tsv>`_, and copy it (Edit->Copy or cntl-c).

2. Paste this into `Tab Separated Values` field on the `Data Upload
   Page </dismod/data/upload>`_.  (Alternatively, you can save this
   file on your computer and upload it as a file, which is good for files
   too large to cut-and-paste.)

3. Click `Create` to load this data into DisMod.  This will take you
   to the model summary page.  The specks of color to the left of the
   table are a tiny summary of all the data you just loaded, in a format that
   Tufte might call the small-multiples sparkplot.  Try clicking on a cell in
   the sparkplot to see a detailed view of some of the
   data you just loaded.

4. Click on `Adjust Priors` to set priors on the data.  You should
   always start by adjusting the priors, but not changing anything.  This
   ensures you know what the DisMod defaults lead to.

5. Click `Apply` to keep the default priors (Note: this only works perfectly on Windows Firefox, for other browsers the java applet does not always redirect to the new page.  To work around this, click on the `View/Modify` link and find your newly cloned model).

6. Click `Calculate covariates for model data`
   to load the age weights for the data set.  Note, this may take 1-2
   minutes, because it is interpolating population data for country-year
   specific population structure.


7. Click `Estimate empirical priors` to fit global pooled data for
   each parameter type separately.  Select `Yes` to generate posterior
   estimates when the empirical prior phase is complete, and select the
   region `Asia, Southeast` to generate posterior estimates for.

8. Be patient; it should take about 1 minute for DisMod to fit the
   empirical priors for this model.  Then it will take an additional
   10 minutes for DisMod to fit the region-sex-year specific posteriors.

9. When DisMod pops up a dialog box to alert you that the posterior
   estimation has completed, click `Summary`, and select a spark plot
   with posterior information (denoted by a solid blue curve) to see how
   the fit looks.  It will look ok, but there is room for improvement.

10. To improve the fit, select `Adjust Priors`, and in the `Adjust Priors` control panel:

   a. Set a level bound on relative risk of 3 for ages before 99
   b. Click `Apply`.

11. Go to step 6, and repeat this process, making just one or two
    changes at a time until you have a satisfactory fit.
