#!/usr/bin/env Rscript
suppressMessages(library("gap")) # for qqplot

# Keep uncommented the set to plot 
#     "HLPHMP" has 10 seed, each with 1000 pairs for a total of 10,000 - to test p-val = 1e-04
#     The others have 2 seeds, each with 500 pairs for a total of 1000
opt = list(id = "HLPHMP", 
           seeds = "17.18.19.20.21.22.23.24.25.26") 
## opt = list(id = "HLPBLL", 
##            seeds = "15.16") 
## opt = list(id = "NKiyTS",
##            seeds = "15.16")
## opt = list(id = "BIrHMP",
##            seeds = "15.16") 
## opt = list(id = "BAPBCr",
##            seeds = "15.16") 

# NB: because of software reasons, to test corr_Ad2s1 ≠ 0 we had to swap traits
#     here trait 2 = trait affected by IGE, trait 1 = trait affected by DGE (proxy)
corr0 = "corr_Ad1s2" 
altfile = "~/PRJs/CFW/output/VDresults/noBatch_pruned_dosages_include_DGE_IGE_cageEffect_corrAd1s2_alt_all.txt"
# Get options 
id = opt$id 
seeds = as.numeric(na.omit(unlist(strsplit(opt$seeds, '[.]'))) )
print(seeds)

# Define input and output depending on seed and id of each set
sim_files = paste0("~/nf_PRJs/nf-CoreQuantGen/simulations/output/simCFW/VD/cage_623/set",id,"/DG1_IG2/0.0/",seeds,"/bivariate/mockphenos/pruned_dosages_DGE_IGE_cageEffect_corr_Ad1s2_all_est.txt")
outfile = paste0("./plot/CFW/Sfig2.qqplot_", id, ".pdf")

# Read realdata to get the names of phenotypes and the p-value observed in real data
VCs.mat = read.table(altfile, header = T, sep="\t")
phenopair = VCs.mat[VCs.mat$pairsID == id, "taskID"] 
pval_obs = VCs.mat[VCs.mat$taskID == phenopair, "pv_chi2dof1"] 

# Define function to parse results
prepare_res <- function(res, nocol = c('sample_size','sample_size_all','covariates_names', 'conv', 'LML')){
  # assigning taskID = name phenoytpe1_phenotype2
  #     TODO: check how it works with univariate, imagine it would be something like trait1_None, trait2_None...
  #     if want something different can do 
  #     if ("trait2" %in% nocol) {res[,"taskID"] = res[,"trait1"] # or res[,"taskID"] = paste0(res[,"trait1"], "univariate")}
  order_cols = c("taskID",colnames(res))
  res[,"taskID"] = paste(res[,"trait1"], res[,"trait2"], sep='_') # NB: like this also have the same name that is saved as _est.txt _STE.txt output files
  res = res[,order_cols]
  
  # filtering out unwanted columns
  res = res[,! colnames(res) %in% nocol]
  # removing columns with all -999 i.e. NAs
  NAs <-  apply(res, 2, FUN = function(res) all(res == -999))
  res = res[!NAs]
  res[res==-999] = NA # see if want to try this too
  return(res)
}
colest = c('trait1', 'trait2', 'sample_size1', 'sample_size1_cm', 'sample_size2', 'sample_size2_cm', #6
           'union_focal', 'inter_focal', 'union_cm', 'inter_cm', #4
           'covariates_names', 'conv', 'LML', #3
           'prop_Ad1', 'prop_Ad2','prop_As1', 'prop_As2', #4
           'corr_Ad1d2', 'corr_Ad1s1', 'corr_Ad1s2', 'corr_Ad2s1','corr_Ad2s2', 'corr_As1s2', #6
           'prop_Ed1', 'prop_Ed2','prop_Es1', 'prop_Es2', #4
           'corr_Ed1d2', 'corr_Ed1s1', 'corr_Ed1s2', 'corr_Ed2s1', 'corr_Ed2s2', 'corr_Es1s2', #6
           'prop_Dm1', 'prop_Dm2', 'corr_Dm1Dm2', #3
           'prop_C1',  'prop_C2', 'corr_C1C2', #3
           'tot_genVar1', 'tot_genVar2', #2
           'total_var1', 'total_var2') #2

# Getting results from all simulations analysed
sims <- lapply(sim_files, function(sim_file) {
  sim = read.csv(file = sim_file, sep = "\t", header = F)
  if(length(colest) != ncol(sim)){
    stop("Something wrong with number of columns and colnames")
  }
  colnames(sim) = colest
  sim = prepare_res(sim, nocol = c(""))
  return(sim)
})
sims = do.call(rbind, sims)
#return(sim)
dim(sims)

row_constrained = which(sims[,corr0] == 0)
names_constrained = sims[row_constrained,"taskID"]
names_full = sims[-row_constrained,"taskID"]
if(!all(names_constrained %in% names_full)) stop("don't have a full model for all corr constrained")
# keeping only the ones that have the constrain (and ordering as row1:full, row2:constrained)
sim_null = sims[row_constrained,]
sim_alt = sims[-row_constrained,]

rownames(sims)[row_constrained] = paste0("null_", sims[row_constrained, "trait1"])
rownames(sims)[-row_constrained] = paste0("alt_", sims[-row_constrained, "trait1"])


# Q-Q plot
cat("saving file to: ", outfile, "\n")
pdf(outfile, h=6, w=5)
par(mar = c(5.1, 5.1, 2.1, 1.1))
pval_plot = pchisq(2*(-sim_alt$LML+sim_null$LML),df=1, lower.tail = FALSE)

r = qqunif(pval_plot, ci=T,las=1, main = phenopair, cex.main = 0.9, 
           sub = paste0('P values null ', length(pval_plot),' simulations'), pch = 18, cex=1.5,
           cex.axis = 1.25, cex.lab=1.4, cex.sub=0.9, col = "#4DBBD5FF", lcol = "#E64B35FF")
qqunif(pval_plot, ci=T,las=1, main = paste0("p-val = ", round(pval_obs, 5)), cex.main = 0.9, 
       sub = paste0('P values null ', length(pval_plot),' simulations'), pch = 18, cex=1.5,
       cex.axis = 1.25, cex.lab=1.4, cex.sub=0.9, cex.main = 1.6, col = "#4DBBD5FF", lcol = "#E64B35FF", 
       #ylim=c(0,5.1), xlim = c(0,5.1))
       ylim=range(r$x,r$y), xlim = range(r$x,r$y))
dev.off()
