This file contains a short description of our implementation.

The most important files that can be run to produce the main results are: "Final_resuts.R" and "Train_test_comparison.R".

There are 4 folders:

1) Analysis
It contains all the files that produce the results.
- "Data_Analysis.R" runs our model with final parameters on the entire data and creates the excel table with the final coefficients. (The table and plots are already in folder Results, but can be created again from scratch by running this file.) Running this file might take around 3 hours. THIS THE MAIN FILE FOR THE PAPER. IT HAS THE MOST IMPORTANT RESULTS.  
- "Grid_search_local.R" contains an example on how grid search can be performed. 
- "Train_test_comparison.R" contains the prediction on test set with best Random forest model, Lasso model, and SH Lasso model. The RF and Lasso models are trained for best coefficients, while the best SH Lasso coefficients are loaded .(SH Lasso can also be trained from scratch, but this might take 2 hours.) The parameter tuning was done with cross validation on train set. 
- "Plots_and_matrix_form_comparison.R" contains some plots and the comparison discussed in section 3 in Supplement. It does not contain the main results of the paper.

2) coefficients
It contains rds files with the coefficients obtained by using our model for given penalizations on all data/train-test split of our data. The names of the files contain the regularization term (lambda) for the main effects, two-way interactions and three-way interactions in this order. 

3) libs
It contains the files with all the functions used to build our model and generate the results.
-"Helpers.R" contains functions used to plot and evaluate the predictions.
-"Combinatorics.R" contains the functions that make possible to identify specific interactions and main effects and use them in our model.
-"Functions_for_updates.R" contains the functions that are used to minimize the penalized negative log-likelihood.
-"Updates.R" contains the algorithms that update the coefficients of our model. It uses "Combinatorics.R" and "functions_for_updates.R". 
-"SHIM_4way_CB.R" contains our continuous Bernoulli hierarchical model. It uses "Combinatorics.R", "functions_for_updates.R", "Updates.R" and "Helpers.R"
-"Read_BHA_data.R" contains a function that loads the Buchwald-Hartwig amination data in the suitable form for our model.
-"Yield_table_full_data.R" contains a function that creates an excel table with results for a given set of coefficients. It will be stored in "Results" folder.

4) Results
-It contains the table of coefficients obtained by using our hierarchical model on the entire data.
-It also contains the plots for predicted vs observed for Lasso, Random Forest, SH Lasso on test data and SH Lasso on entire data.

Information about running the code:

-Running "Final_resuts.R" and "Train_test_comparison.R" will produce the main results from the paper and the supplement of the paper.

-We will report the running times for our model on an Intel(R) Core(TM) Ultra 7 155U 1.70 GHz. Training our hierarchical model from scratch on the entire data takes between 1 and 15 hours depending on the amount of regularization considering tolerance level for convergence between 5*10^(-2) and *10^(-3). The higher the regularization, the faster is the model due to the hierarchical constraints. Training the model on entire data for the final amount of penalization selected by AIC takes 3 hours.     
