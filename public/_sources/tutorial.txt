=====================
Quick-start Tutorials
=====================

Best-case scenario:  Lots of consistient, low-noise data
--------------------------------------------------------

`This Data Table <asthma_data.html>`_ is unrealistically good data, but it provides a
starting point for learning how to use DisMod.

1. Select all the rows and columns from the `asthma data table <asthma_data.html>`_, and copy it (Edit->Copy or cntl-c)
2. Paste this into `Tab Separated Values` field on the `Import Rates Page </rate/>`_
3. Click 'Create' to load this data into DisMod.
4. Try fitting one of these age-specific rate functions individually.
   For example, if the id number of the first plot is 957, as shown in this screenshot,

   .. image:: images/tutorial_1a.png

   type this on the command line::

    scripts/fit_rf 957
5. This is pretty good, but it seems to have a strange behavior at the oldest ages:

   .. image:: images/tutorial_1b.png

   To add prior beliefs that the incidence rate is a smooth function
   of age and is decreasing at older ages, click 'clone' next to
   957. Then in the 'Priors' textbox, add the lines::

    smooth 10.0
    decreasing 80 100

   Now, the results of fitting the new rate function (with
   ``scripts/fit_rf <new_rate_id>``) incorporate these prior beliefs:

   .. image:: images/tutorial_1c.png

   It might be a good idea to add similar priors to the remission rate.

6. To fit all the data together (and generate incidence and duration
   estimates as well), do the following (but replace the appropriate rate id
   numbers based on your situation)::

    scripts/fit_dm 959 960 961

   .. image:: images/tutorial_1d.png

7. To export the estimates for use in a spreadsheet or another
   program, direct your browser to the rate function with .csv instead
   of .html in the url.


Another example, with more, noiser, less consistent data
--------------------------------------------------------

Now try the same thing with `example data table <example_data.html>`_,
which is some synthetically generated data that is similar to the
asthma data from the previous section.
