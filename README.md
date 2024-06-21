This repository contains code and data used in the master's thesis "Inference and Prediction in Stochastic Dynamical Systems" by Axel Ming Szetu.
The code was developed in R version 4.3.2.

Several of the notebooks use functions defined in "AMOCestimaton.Rmd" from the following repository:
https://erda.ku.dk/archives/cb78329f209d8ff2b4dd810abe4780ae/published-archive.html
For this reason, "AMOCestimaton.Rmd" should be run in the envieronment before running any of the notebooks in this repository.

The code for associated with Chapter 3 is found "ROC curves for Variance and Autocorrelation.Rmd".
It uses the data files "X.H0.Rdata", "X.H1.Rdata" and "EWS.matrices.new.windows.Rdata".
While the default setting is to load these data sets from storage, the code to generate similar datasets is included and can be run.
The simulations were not seeded.

The code associated with Chapter 4 is found in "Misspecification of Drift.Rmd"

I apologize for the mess in the folder.
