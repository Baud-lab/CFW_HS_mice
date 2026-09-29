# Data results from bivariate analysis to understand the mechanism of indirect genetic effects

In each population folder (`CFW` and `HSmice`) there are:

+ `dataset`: 
  + h5 file: storing all information for VD analysis (genotypes, phenotypes, 
  covariates, cages) 
  + file: with information on macrophenotype and category for all phenotypes
  
+ `permuations`: 
  + results from VD analysis of the most significant pair of phenotypes while scrambling 
  cage assimgnments to model IGE

+ `simulations`:
  + `params`: values used to simulate pairs of phenotypes (one file per set)
  + results from VD analysis of simulations `bivariate` and `sexvariate`

+ `VD`:
  + results from VD analysis on real data 



