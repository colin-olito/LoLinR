# local-lin PROJECT SUMMARY AND FILE DESCRIPTIONS: 

##############
#  SUMMARY
This is a project with the aim of developing a local linear regression function for use with metabolic time-series data. The primary motivation for this function was to design a statistically robust and repeatable method for estimating resting metabolic rate from O2 consumption time series data. To do this, I modified a basic LOESS regression technique. The function starts by performing local linear regressions on increasingly large subsets of the full input data set. If we are interested in just the slope of these regressions, then the current function is functional.  

My main worry is that this function alone is not substantial enough to warrant a methods paper, even if it would be generally useful for physiologists and ecologists interested in metabolism, and with access to the equipment necessary to generate these kinds of time-series.  However, I do think the general approach is robust and flexible... so I really wonder if it wouldn't be relatively easy to generalize the technique to calculate any number of potentially useful metrics for these kinds of data.  I know Craig mentioned that it is sometimes necessary to determine what proportion of a data set should be used for estimating something from flow-through respirometry studies... perhaps it wouldn't take much to expand the code base here to provide a "toolkit" of particularly useful metabolic time-series R functions. If we did, then I think it would make more sense to turn the collection of functions into an R package - something Diego has offered his help with. If we did this, I think it would be much easier to justify the usefulness of our contribution to a bigger methods journal like Methods in Ecology and Evolution...  

################
#  Package outline

- Reading Data Files
   - Need to think of file types (i.e. - Firesting output .txt's.


- Plotting 
   - FindLocLin: 'plots' option creates panel of best 25 subsets.
   - PlotBest(): Plots resids from 1st ranked subset, or user defined.
   - distribution of L.

TO DO:
   - Defensive programming for obvious errors
   - Option for user specified window sizes?
   - Diagnostics for combined metric. e.g. L ~ window size.
   - 'Export' everything



################
#  Testing

- Obvious tests
   - math problems
   - defensive errors
   - Alignment b/w FindLocLin() & PlotBest() - Make sure these functions return the exact same regression results.


- Corner cases
   - When to expect errors
   - 




################
#  File Descriptions

custom_loess_devel.R -- scratch code I wrote while initially developing the FindLocLin function

FindLocLin_functions.R -- Clean, working code with the necessary functions, and explanations

testplots.pdf -- Scatterplots of the best 25 local regressions, giving an idea of what the function does

residplots.pdf -- residual plots of the 'best' linear regression

taxes.csv -- test data from an old stats course I took... useful for test-driving the function because it's ugly and non-linear

TestO2data.csv -- Test O2 consumption data for the sea urchin H. erythrogramm from my empirical work last year.



