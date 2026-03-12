est_wME = "~/PRJs/HSmice/output/VDresults/data_bcNcovariates_Andres_kinship_None_DGE_IGE_cageEffect_maternalEffect_corr_Ad1s2_alt_all.txt"
est_noME = "~/PRJs/HSmice/output/VDresults/data_bcNcovariates_Andres_kinship_None_DGE_IGE_cageEffect_corr_Ad1s2_alt_all.txt"
outpdf = "./plot/HSmice/SFig5.corr_wtME_vs_noME.pdf"

wtME_alt = read.delim(est_wME)
rownames(wtME_alt) = wtME_alt$taskID

noME_alt = read.delim(est_noME)
rownames(noME_alt) = noME_alt$taskID

# Now filter for wtME
noME_alt = noME_alt[rownames(wtME_alt),]

pdf(outpdf, h=7, w=7)
par(mar = c(5.1, 5.1, 4.1, 1.1))
sapply(grep("corr_Ad2s1",colnames(noME_alt), value=T), 
       function(p){cor = round(cor.test(wtME_alt[,p], noME_alt[,p])$estimate, 2)
       cor_sig = ifelse(cor.test(wtME_alt[,p], noME_alt[,p])$p.value < 0.05, "<", ">")
       plot(x = wtME_alt[,p], y = noME_alt[,p], 
            main = paste0("correlation = ", cor, " (p-value ",cor_sig, " 0.05)"), cex.main = 0.8,
            ylim = c(-1,1), xlim=c(-1,1),
            xlab = paste0(p, " (DGE+IGE+CE+ME)"),  ylab = paste0(p," (DGE+IGE+CE)"), 
            cex.lab = 1.4, cex.axis=1.2);
       abline(a=0, b=1, lty=2, col="grey")})
dev.off()
