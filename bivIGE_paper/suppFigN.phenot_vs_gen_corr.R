sourcefun= "/users/abaud/htonnele/git/me/core_VD/code/Rsrc/functions/" # TODO: might have to change this path relative to where code is
source(file.path(sourcefun, "select_col_VCsmat.R")) # `select_col`: function to select columns of a specific Variance component, for all phenotypes 
sourcemod= "/users/abaud/htonnele/git/me/core_VD/code/Rsrc/afterVD/real/modules/"

# Corr_Ads estimated in univar vs corr_Ad1s1 or corr_Ad2s2 estimated in bivar
pop = "CFW"
bi_est_dir = "/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/realdata/output/CFW/VD/bivariate/noBatch/pruned_dosages_include_DGE_IGE_cageEffect/"
macro_file = "/users/abaud/htonnele/PRJs/CFW/output/dataset/macropheno_CFW.csv"
h5file = "/users/abaud/htonnele/PRJs/CFW/output/dataset/CFWmice_HT.h5"
phenoV = "noBatch"
outplot = "./plot/CFW/SFigN.phenot_vs_gen_corr.pdf"
selected_pairs = list(DGE = c("WH.Ears_Area", "Bioch.Chloride_BWcorr", "Bioch.Glycerol_BWcorr", "Bioch.LDL_BWcorr", "Haem.MPV", "Haem.Large_PLT", 
                              "Bioch.Glucose_BWcorr", "FACS.CD3posCD4posCD44pos", "FACS.CD3posCD8posCD44pos", "BMC.osteoporosis", "Haem.abs_neuts", 
                              "Hypoxia.f_SHR_BWcorr", "Hypoxia.f_Baseline_BWcorr", "Hypoxia.f_Undershoot_BWcorr"),
                      IGE = c("Adrenals.Adrenals_g_BWcorr", "Bioch.LDL_BWcorr", "Haem.MPV", "Haem.Large_PLT", "PST.Immobility.Last4min", 
                              "PST.Immobility.First2min",  "Bioch.Amylase_BWcorr", "Haem.EOS_percent", "FACS.CD3posCD44negCD4CD8Ratio", 
                              "FACS.CD3posCD4CD8Ratio", "FACS.CD3posCD8pos"))

## pop = "HSmice"
## bi_est_dir = "/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/realdata/output/HSmice/VD/bivariate/data_bcNcovariates/Andres_kinship_None_DGE_IGE_cageEffect/"
## macro_file = "/users/abaud/htonnele/PRJs/HSmice/output/dataset/macropheno_HSmice.csv"
## h5file = "/users/abaud/htonnele/PRJs/HSmice/output/dataset/HSmice_wPleth_noEpoch_HT_v2.h5"
## phenoV = "data_bcNcovariates"
## outplot = "./plot/HSmice/SFigN.phenot_vs_gen_corr.pdf"
## selected_pairs = list(DGE = c("Imm.BCell.size", "Glucose_0", "Biochem.HDL", "CD8Count", "CD4Count", "Haem.LYMabs", "Haem.WBC", "FN.preWeight", "FN.postWeight", "Start.Weight", "Weight.GrowthSlope", 
##                               "Imm.CD4inCD3XGeoMean", "AdrenalMeanWeight", "Imm.CD8inCD3YGeoMean", "Insulin.75", "Pleth.base.EnhancedPause", "Biochem.Calcium", "Glucose_15"), 
##                       IGE = c("Imm.CD4inCD3XGeoMean", "Imm.BCell.size", "Imm.CD4.size", "Imm.CD8.size", "Biochem.ALT", "Glucose_0", "Biochem.Calcium", "Biochem.Chloride", 
##                               "Biochem.Sodium", "Biochem.Creatinine", "Biochem.Triglycerides", "Biochem.HDL", "Glucose_15", "CD4Count", "Haem.LYMabs", "Haem.WBC", "Pleth.base.TidalVolume", 
##                               "Weight.GrowthSlope", "Start.Weight", "Pleth.base.EnhancedPause", "FN.preWeight", 
##                               "FN.postWeight", "Obesity.BodyLength", "Cue.Mean.Freeze.Corrected.During"))


# Starting here analysis
files = list.files(bi_est_dir, pattern = "_estNste.Rdata")
source(file.path(sourcemod, "6.1.parseData.R")) 

# getting names of trait1 and trait2 from VC matrix - going to need this to order and filter matrices
trait1 = gsub(".trait1", "", grep("trait1", colnames(VCs.mat), value = T))
trait2 = rownames(VCs.mat)
print("dim VCs.mat: ")
print(dim(VCs.mat))
trait1N2 = unique(c(trait1, trait2)[order(factor(c(trait1, trait2), levels=macropheno))])

corr_Ad2s1 = select_col(VCs.mat, "corr_Ad2s1")
colnames(corr_Ad2s1) = sub(".corr_Ad2s1", "", colnames(corr_Ad2s1))


## get phenotypic correlation
suppressMessages(library("Hmisc"))
if(is.null(h5file)) stop("Need h5file to get phenotypic measurements")
if(is.null(phenoV)) stop("Need pheno version to get phenotypic measurements")
# 2. cor_mat and p_mat FOR PHENOTYPIC CORS -----------
get_pheno = \(){
  ###### pheno from h5
  pheno = h5read(h5file, paste0("/phenotypes/",phenoV,"/matrix")); 
  colnames(pheno) = h5read(h5file, paste0("/phenotypes/",phenoV,"/col_header/phenotype_ID"));
  rownames(pheno) = h5read(h5file, paste0("/phenotypes/",phenoV,"/row_header/sample_ID"));
  pheno[which(pheno == -999)] = NA
  pheno = as.data.frame(pheno[!apply(pheno, 1, function(x) all(is.na(x))),])
  return(pheno)
}

pheno = get_pheno()

### in this way they are ALL PEARSON
###    Get correlation between phenotypes 
## res = rcorr(as.matrix(pheno), type = "pearson")
## diag(res$P) = 0 # setting the diagonal p-value to zero
## 
## cor_mat = res$r[trait1N2, trait1N2]
## #    Get Pvalue of correlation between phenotypes
## p_nom = res$P[trait1N2, trait1N2]


### Function getting the correlation between two phenotypes, either pearson or spearman, 
choose_method = function(x, y) {
  # Test for normality
  p1 = shapiro.test(x)$p.value # if > 0.05, normally distributed
  p2 = shapiro.test(y)$p.value # if > 0.05, normally distributed
  
  # Test for linearity using Pearson vs. Spearman
  pearson = cor.test(x, y, method = "pearson", use = "pairwise.complete.obs") # selecting only non-NAs in both 
  spearman = cor.test(x, y, method = "spearman", use = "pairwise.complete.obs", exact = F) # selecting only non-NAs in both 
  
  # Compare Spearman and Pearson results # Checking for linearity
  linear_check = abs(pearson$estimate - spearman$estimate) < 0.1  # Threshold for linearity - heuristic threshold
  
  # Choose method
  if (p1 > 0.05 && p2 > 0.05 && linear_check) {
    # both phenotypes are normally distributed and they have a linear relathionship
    # hence using pearson
    return(list(method = "pearson", corr = pearson$estimate, pval = pearson$p.value) )
  } else {
    # using spearman when one of the above conditions is not met
    return(list(method = "spearman", corr = spearman$estimate, pval = spearman$p.value) )
  }
}

# Initialize correlation and method matrices
n = ncol(pheno)
corr_mat = matrix(NA, nrow = n, ncol = n)
met_mat = matrix("", nrow = n, ncol = n)
colnames(corr_mat) = colnames(pheno)
rownames(corr_mat) = colnames(pheno)
colnames(met_mat) = colnames(pheno)
rownames(met_mat) = colnames(pheno)

# Pairwise correlation calculation
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    if (!sum(complete.cases(pheno[, c(i,j)])) >= 2) next
    x = pheno[, i]
    y = pheno[, j]
    
    # Choose the method
    corr_xy = choose_method(x, y)
    
    # Calculate correlation
    #correlation = cor(x, y, method = method, use = "pairwise.complete.obs")
    
    # Store results
    corr_mat[i, j] = corr_mat[j, i] = corr_xy$corr
    met_mat[i, j] = met_mat[j, i] = corr_xy$method
  }
}
diag(corr_mat) = 1
pheno_mat = corr_mat

pheno_mat = pheno_mat[rownames(corr_Ad2s1), colnames(corr_Ad2s1)]
dim(pheno_mat)

stopifnot(all.equal(dimnames(pheno_mat), dimnames(corr_Ad2s1)))

coolors = matrix("black", nrow = nrow(corr_Ad2s1), ncol= ncol(corr_Ad2s1), dimnames = dimnames(corr_Ad2s1))
coolors[selected_pairs$DGE, selected_pairs$IGE] = "darkorange"
pch = matrix(1, nrow = nrow(corr_Ad2s1), ncol= ncol(corr_Ad2s1), dimnames = dimnames(corr_Ad2s1))
pch[selected_pairs$DGE, selected_pairs$IGE] = 16

pdf(outplot)
par(mar=c(5.1, 5.1, 4.1, 1.1))
plot(as.matrix(corr_Ad2s1), pheno_mat, xlab = "IGE-DGE correlation", ylab = "phenotypic correlation", cex.axis = 1.2, cex.lab = 1.4, 
     col = coolors, pch = pch, xlim = c(-1,1), ylim = c(-1,1))
abline(a = 0, b = 1, col = "red", lty = 2)
mtext(line = 1, paste0("cor = ", round(cor.test(c(pheno_mat), c(as.matrix(corr_Ad2s1)))$estimate, 2), 
                       "; p-val < 2.2e-16"))
#points(as.matrix(corr_Ad2s1[selected_pairs$DGE, selected_pairs$IGE]), pheno_mat[selected_pairs$DGE, selected_pairs$IGE], 
#     col = coolors[selected_pairs$DGE, selected_pairs$IGE], pch = pch[selected_pairs$DGE, selected_pairs$IGE])
dev.off()
# Checking distributions, NB: phenotypic corr has some values = 1 - the self phenotype
hist(c(pheno_mat))
hist(c(as.matrix(corr_Ad2s1)))


