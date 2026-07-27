### CFW ####
sexres = "/users/abaud/htonnele/PRJs/CFW/output/nf-output/VD/sexvariate/noBatch/pruned_dosages_include_DGE_IGE_IEE_cageEffect_all_estNste.Rdata"
#unires = "/users/abaud/htonnele/PRJs/CFW/output/VDreal_2310/univariate/noBatch_500/pruned_dosages_include_DGE_IGE_cageEffect_estNste.Rdata"
macro = "/users/abaud/htonnele/PRJs/CFW/output/dataset/macropheno_CFW.csv"
selected_pairs = list(DGE = c("WH.Ears_Area", "Bioch.Chloride_BWcorr", "Bioch.Glycerol_BWcorr", "Bioch.LDL_BWcorr", "Haem.MPV", "Haem.Large_PLT", 
                              "Bioch.Glucose_BWcorr", "FACS.CD3posCD4posCD44pos", "FACS.CD3posCD8posCD44pos", "BMC.osteoporosis", "Haem.abs_neuts", 
                              "Hypoxia.f_SHR_BWcorr", "Hypoxia.f_Baseline_BWcorr", "Hypoxia.f_Undershoot_BWcorr"),
                      IGE = c("Adrenals.Adrenals_g_BWcorr", "Bioch.LDL_BWcorr", "Haem.MPV", "Haem.Large_PLT", "PST.Immobility.Last4min", 
                              "PST.Immobility.First2min",  "Bioch.Amylase_BWcorr", "Haem.EOS_percent", "FACS.CD3posCD44negCD4CD8Ratio", 
                              "FACS.CD3posCD4CD8Ratio", "FACS.CD3posCD8pos"))
outpdf = "./plot/CFW/SFigN.dotplot_sexest_igetrait.pdf";h=3.5;w=6.5

### HSmice ####
sexres = "/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/realdata/output/HSmice/VD/sexvariate/data_bcNcovariates/Andres_kinship_DGE_IGE_IEE_cageEffect_corr_As1s2_one_estNste.Rdata"
#unires = "/users/abaud/htonnele/PRJs/HSmice/output/VDreal_2311/univariate/data_bcNcovariates_500/Andres_kinship_DGE_IGE_cageEffect_estNste.Rdata"
macro = "/users/abaud/htonnele/PRJs/HSmice/output/dataset/macropheno_HSmice.csv"
selected_pairs = list(DGE = c("Imm.BCell.size", "Glucose_0", "Biochem.HDL", "CD8Count", "CD4Count", "Haem.LYMabs", "Haem.WBC", "FN.preWeight", "FN.postWeight", "Start.Weight", "Weight.GrowthSlope", 
                              "Imm.CD4inCD3XGeoMean", "AdrenalMeanWeight", "Imm.CD8inCD3YGeoMean", "Insulin.75", "Pleth.base.EnhancedPause", "Biochem.Calcium", "Glucose_15"), 
                      IGE = c("Imm.CD4inCD3XGeoMean", "Imm.BCell.size", "Imm.CD4.size", "Imm.CD8.size", "Biochem.ALT", "Glucose_0", "Biochem.Calcium", "Biochem.Chloride", 
                              "Biochem.Sodium", "Biochem.Creatinine", "Biochem.Triglycerides", "Biochem.HDL", "Glucose_15", "CD4Count", "Haem.LYMabs", "Haem.WBC", "Pleth.base.TidalVolume", 
                              "Weight.GrowthSlope", "Start.Weight", "Pleth.base.EnhancedPause", "FN.preWeight", 
                              "FN.postWeight", "Obesity.BodyLength", "Cue.Mean.Freeze.Corrected.During"))
outpdf = "./plot/HSmice/SFigN.dotplot_sexest_igetrait.pdf";h=6;w=7


# load results
load(sexres)
VCs = res$VCs
rownames(VCs) = gsub("_sex1", "",VCs[,"trait1"])

## # Select trait_ige
## load(unires)
## uni_VCs = uni_VCs[gsub("_sex1", "", VCs[,"trait1"]),]
## uni_VCs = uni_VCs[order(uni_VCs[,"prop_As1"], decreasing = F),]
## trait_ige = uni_VCs[uni_VCs[,"prop_As1"] > 0.05,"trait1"]
## trait_ige = trait_ige[trait_ige != "Adrenals.Adrenals_g"]

# Select trait_ige for only in cluster
trait_ige = selected_pairs$IGE
sel = VCs[trait_ige[length(trait_ige):1],]
sel[abs(sel[, "corr_As1s2"]) < 0.5, c("trait1", "corr_As1s2")]
#sel = sel[order(sel$corr_As1s2, decreasing = T), ]

# reaading macrophenotypes to set colours
dict_ord = c("immunology", "blood", "glucose", "steroid","physiology", "anthropometric", "brain", "other", "behaviour") # behaviours last
ormacro = read.csv(macro, header = T)
ormacro = ormacro[order(factor(ormacro$category, levels=dict_ord)), ]
table(ormacro[,"category"])
macropheno = setNames( ormacro$phenotype_ID, ormacro$macrophenotype)
cat = setNames( ormacro$phenotype_ID, ormacro$category)

### specific CFW
if (length(grep("CFW", macro)) != 0){
  replacements <- c("Adrenals" = "AdrenalWeight", 
                    "Bioch" = "Biochemistry", 
                    "WH" = "EarPunch",
                    "Haem" = "Haematology", 
                    "FC" = "Cue")
  names(macropheno) = sapply(names(macropheno), 
                             function(x, y) ifelse(x %in% names(y), return(y[x]), x), 
                             replacements, USE.NAMES = F)
}

coolors = c("immunology"     = "#53585F", 
            "blood"          = "#5D6D84", 
            "glucose"        = "#8DA0CB",
            "steroid"        = "#9EC2C9",
            "physiology"     = "#8ABFAA",
            "anthropometric" = "#A4CBA6",
            "bones"          = "#88A77D",
            "brain"          = "#5C8A72", 
            "other"          = "#004242",
            "behaviour"      = "#B894B1")

points_vc <- function(par_name, sel, cex.pheno=0.8){
  ste_dict = gsub("corr_", "STE_", gsub("prop_", "STE_", grep("corr_A|prop_A", colnames(sel), value = T)) )
  names(ste_dict) = grep("corr_A|prop_A", colnames(sel), value = T)
  ste_name = ste_dict[par_name]
  
  #sel = sel[order(sel$corr_As1s2, decreasing = T), ]
  est = sel[,par_name]
  names(est) = rownames(sel)
  ste = sel[,ste_name]
  names(ste) = rownames(sel)
  stopifnot(all(names(est) == names(ste)))
  print(length(est))
  print(length(ste))
  trait_cols = coolors[names(cat[na.omit(match(names(est), cat))])]
  
  Lbord = 15 # NB: it is very impo that the left and right border are the same for plot and histogram
  Rbord = 7 
  xlimi=range(c(est-ste, est+ste)) # need to be the same for dotplot and histogram
  if(length(grep("corr", par_name)) != 0){
    xlimi=range(c(-1, 1)) # need to be the same for dotplot and histogram
  }
  # if want to be a bit larger on the limit, do the ones below
  #ronge=range(c(est-ste, est+ste))
  #sprk = (ronge[2] - ronge[1]) / 10 # need to be the same for dotplot and histogram
  #xlimi= c(ronge[1]- sprk, ronge[2] + sprk)
  
  #### 1. Dot plot ####
  # a. starting with an empty plot - set the axis and labels 
  #    will plot the dot later to have them on top of the arrow lines
  par(mar=c(5, Lbord, 1, Rbord)) # NB: par(mar) = par(mai), the difference is in the unit of measure, 'mar' is in 'lines'; 'mai' is in 'inches'
  
  plot(x = est, y = 1:length(est), 
       xlim=xlimi, xlab=par_name, ylim = range(1:length(est)),
       ylab="", yaxt="n",
       cex.axis=0.8, type="n")
  
  abline(v = 0, lty=5, col="grey90", lwd=1.5) 
  
  # b. adding horizontal labels with pheno
  abline(h = 1:length(est), col = "grey90", lty = 3) 
  axis(at = 1:length(est), side = 2, las =1, cex.axis=cex.pheno, lwd=0, lwd.ticks = 1, 
       labels = names(est))
  
  # c. defining color and dot type depending on values
  #pch_dot = dot_col = cex_dot = vector(length = length(est)); 
  dot_col = trait_cols #"black"
  pch_dot = 16
  cex_dot = 1
  
  # d. adding the +/- STE, as arrows with center the estimate
  arrows(est-ste, 1:length(est), est+ste, 1:length(est), 
         length=0.02, angle=90, code=3, col=adjustcolor(dot_col, 0.8), lty=1) # code=3: head at both ends of the arrow
  # e. adding the points! 
  points(est, 1:length(est), pch=pch_dot, col=dot_col, cex = cex_dot) 
  
  #mtext("pheno", side=2, line=13.5)
  
  legend(x = par()$usr[2], y = par()$usr[4], lty = rep(1, length(unique(trait_cols))), lwd = rep(5, length(unique(trait_cols))), 
         col = unique(trait_cols), legend = unique(names(trait_cols)), 
         xpd = T, bty = "n", cex = 0.8, x.intersp = 0.5, seg.len = 0.5)
  
}


pdf(outpdf, h=h, w=w)
# no histogram
points_vc("corr_As1s2", sel)
dev.off()
