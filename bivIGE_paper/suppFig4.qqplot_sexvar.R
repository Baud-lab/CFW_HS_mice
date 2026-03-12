#!/usr/bin/env Rscript
suppressMessages(library("gap")) # for qqplot

real_file = "~/nf_PRJs/nf-CoreQuantGen/realdata/output/CFW/VD/sexvariate/noBatch/pruned_dosages_include_DGE_IGE_IEE_cageEffect_corr_As1s2_one_all_estNste.Rdata"
real_pheno = "Haem.EOS_percent_sex1_Haem.EOS_percent_sex2"
sim_file = "~/nf_PRJs/nf-CoreQuantGen/simulations/output/simCFW/VD/cage_604/setSxHEp/IG1_IG2/1.0/40/sexvariate/mockphenos/pruned_dosages_DGE_IGE_cageEffect_corr_As1s2_one_all_estNste.Rdata"
sim_params = "~/nf_PRJs/nf-CoreQuantGen/simulations/output/simCFW/mockphenos/cage_604/setSxHEp/IG1_IG2/1.0/40/params_sex_V1.0_S40.txt"
outpdf = "./plot/CFW/SFig4.qqplot_sexvar.pdf"
coroi = "corr_As1s2"

# Loading results from real phenotypes
load(real_file)
real_alt = res$VCs
rm(res)

# Loading results from simulated phenotypes
load(sim_file)
alt = res$VCs
null = res$VCs0
rm(res)

# Loading simulated params 
sim_par = read.table(sim_params)

# Gettnig p-val from real data
pval_obs = real_alt[real_alt$taskID == real_pheno, "pv_chi2dof1"]
pval_obs


# Q-Q plot
pdf(outpdf)
par(mar = c(5.1, 5.1, 2.1, 1.1))
pval_plot = pchisq(2*(-alt$LML+null$LML), df=1, lower.tail = FALSE)

r = qqunif(pval_plot, ci=T,las=1, main = paste0("p-val = ", round(pval_obs, 5)),
           sub = paste0('P values ', length(pval_plot),' simulations'), pch = 18, cex=1.5,
           cex.axis = 1.25, cex.lab=1.4, cex.sub=0.9, col = "#4DBBD5FF", lcol = "#E64B35FF")
# uncomment if want plot that has ylim = xlim
## qqunif(pval_plot, ci=T,las=1, main = paste0("p-val = ", round(pval_obs, 5)), 
##        sub = paste0('P values ', length(pval_plot),' simulations'), pch = 18, cex=1.5,
##        cex.axis = 1.25, cex.lab=1.4, cex.sub=0.9, cex.main = 1, col = "#4DBBD5FF", lcol = "#E64B35FF", 
##        ylim=range(r$x,r$y), xlim = range(r$x,r$y))
dev.off()
