###########################################################
############  Script to plot Supp Fig. 4   ################
# Plot to compare IGE-DGE correlations and their standard #
# errors in the two datasets                              #
###########################################################

#   for HS mice and CFW mice
CFW_file = "./data/CFW/VD/bivariate/noBatch_pruned_dosages_include_DGE_IGE_cageEffect_corrAd1s2_alt_all.txt"
HS_file = "./data/HSmice/VD/bivariate/data_bcNcovariates_Andres_kinship_None_DGE_IGE_cageEffect_corr_Ad1s2_alt_all.txt"

outpdf = "./plot/SFig4.compare_cor_se_HSnCFW.pdf"
#outpdf2 = "./plot/SFig4b.cor_vs_pval_HSnCFW.pdf"

boxplot_vc = function(HS, CFW, main, colHS = "darkorange", colCFW = "#C51B7D", ...){
  par(cex.lab = 1.2)
  boxplot(HS, CFW, 
          xaxt = "n", main = main, 
          col= adjustcolor("white", alpha.f = 0), border=adjustcolor("white", alpha.f = 0), ...)
  axis(1, at = c(1,2), labels = c("HS mice", "CFW mice"))
  points(x = jitter(rep(1, length(HS)), factor = 9), HS, 
         col = colHS, pch = 1)
  points(x = jitter(rep(2, length(CFW)), factor = 4.5), CFW, 
         col = colCFW, pch = 1)
  boxplot(HS, CFW, 
          xaxt = "n", yaxt = "n",
          col = adjustcolor("white", alpha.f = 0.5), 
          boxwex = 0.5, outline = F, bty = "n",
          add=T, ...)
  arrox = c(1.35,1.65)
  arrows(x0 = arrox, x1 = arrox, 
         y0 = c(mean(HS) - sd(HS), mean(CFW) - sd(CFW)), y1 = c(mean(HS) + sd(HS), mean(CFW) + sd(CFW)), 
         col = c(colHS,colCFW), length = 0, lwd = 4) # "grey40"
  points(x = arrox, y = c(mean(HS), mean(CFW)), col = c(colHS,colCFW), cex = 2, pch = 15)
  return(list("HSmice" = c("mean" = mean(HS), "sd" = sd(HS)), 
              "CFW" = c("mean" = mean(CFW), "sd" = sd(CFW)) ))
}

pdf(outpdf, h=6, w = 6)
par(las = 1, cex.axis = 1.2)
set.seed(2)
# bivariate
biCFW = read.csv(CFW_file, header = T, sep = "\t")
biHS = read.csv(HS_file, header = T, sep = "\t")

#boxplot(abs(biHS$corr_Ad2s1), abs(biCFW$corr_Ad2s1), main = "corr_Ad2s1")
cor_dge_ige = boxplot_vc(abs(biHS[,"corr_Ad2s1"]), abs(biCFW[,"corr_Ad2s1"]), "IGE-DGE correlation", ylab = "absolute value")

#boxplot(biHS$STE_Ad2s1, biCFW$STE_Ad2s1, main = "STE_Ad2s1")
se_cor = boxplot_vc(biHS[,"STE_Ad2s1"], biCFW[,"STE_Ad2s1"], "se(corr)", ylab = "value")

dev.off()

### If want to see the relationship of estimate with p-val
## pdf(outpdf2, h=6, w = 6)
## par(las = 1, cex.axis = 1.1, cex.lab = 1.2, cex.main = 1.4)
## sig = biHS$fdr < 0.1
## plot(abs(biHS$corr_Ad2s1[!sig]), -log10(biHS$pv_chi2dof1[!sig]),
##      xlim = range(abs(biHS$corr_Ad2s1)), ylim = range(-log10(biHS$pv_chi2dof1)), 
##      main = "HS mice", xlab = "abs(cor)", ylab = "-log10(p-val)")
## points(x = abs(biHS$corr_Ad2s1[sig]), y = -log10(biHS$pv_chi2dof1[sig]), 
##        col = "red", pch = 16, cex = 1.1)
## 
## sig = biCFW$fdr < 0.1
## plot(abs(biCFW$corr_Ad2s1[!sig]), -log10(biCFW$pv_chi2dof1[!sig]),
##      xlim = range(abs(biCFW$corr_Ad2s1)), ylim = range(-log10(biCFW$pv_chi2dof1)),
##      main = "CFW mice", xlab = "abs(cor)", ylab = "-log10(p-val)")
## points(x = abs(biCFW$corr_Ad2s1[sig]), y = -log10(biCFW$pv_chi2dof1[sig]), 
##        col = "red", pch = 16, cex = 1.1, 
##        main = "CFW mice")
## dev.off()
