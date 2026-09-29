##################################################################
############       Script to plot Supp Fig. 2      ###############
# Plot qqplot to check p-value calibration in bivariate analysis #
# for different sets of simulations based on different real      #
# phenotypic pairs with increasing p-values (done only for CFW)  #
##################################################################

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
altfile = "./data/CFW/VD/bivariate/noBatch_pruned_dosages_include_DGE_IGE_cageEffect_corrAd1s2_alt_all.txt"
# Get options 
id = opt$id 
seeds = as.numeric(na.omit(unlist(strsplit(opt$seeds, '[.]'))) )
print(seeds)

# Define input and output depending on seed and id of each set
sim_files = paste0("./data/CFW/simulations/bivariate/set",id,"_DG2_IG1_0.0_s",seeds,"_estNste.Rdata")
outfile = paste0("./plot/CFW/SFig2.qqplot_", id, ".pdf")

# Read realdata to get the names of phenotypes and the p-value observed in real data
VCs.mat = read.table(altfile, header = T, sep="\t")
phenopair = VCs.mat[VCs.mat$pairsID == id, "taskID"] 
pval_obs = VCs.mat[VCs.mat$taskID == phenopair, "pv_chi2dof1"] 

# Getting results from all simulations analysed
sims = lapply(sim_files, function(sim_file) get(load(sim_file))) 

sim_alt = do.call(rbind, lapply(sims, "[[", "VCs"))
sim_null = do.call(rbind, lapply(sims, "[[", "VCs0"))

stopifnot(all(sim_alt$taskID == sim_null$taskID))

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
