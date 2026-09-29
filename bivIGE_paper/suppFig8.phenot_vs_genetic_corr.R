###########################################################
############  Script to plot Supp Fig. 8   ################
# Plot phenotypic correlation averaged over cage mates vs #
# IGE-DGE correlation                                     #
# Plot phenotypic correlation vs DGE-DGE correlation      #
###########################################################

# Libraries for phenotypic correlation
suppressMessages(library("Hmisc"))
suppressMessages(library("rhdf5"))

sourcefun = "./Rfun/"
source(file.path(sourcefun, "select_col_VCsmat.R")) # `select_col`: function to select columns of a specific Variance component, for all phenotypes 

# To plot CFW - comment lines for HS mice
pop = "CFW mice"
bi_est_dir = "./data/CFW/VD/bivariate/noBatch_pruned_dosages_include_DGE_IGE_cageEffect/"
macro_file = "./data/CFW/dataset/macropheno_CFW.csv"
h5file = "./data/CFW/dataset/CFWmice.h5"
phenoV = "noBatch"
cageV = "all623"
outplot = "./plot/CFW/SFig8.phenot_vs_gen_corr.pdf"
pthr = 0.00031 # threshold of nominal p-value based on FDR < 0.1

# To plot HSmice - comment lines for CFW
## pop = "HS mice"
## bi_est_dir = "./data/HSmice/VD/bivariate/data_bcNcovariates_Andres_kinship_None_DGE_IGE_cageEffect/"
## macro_file = "./data/HSmice/dataset/macropheno_HSmice.csv"
## h5file = "./data/HSmice/dataset/HSmice.h5"
## phenoV = "data_bcNcovariates"
## cageV = "all"
## outplot = "./plot/HSmice/SFig8.phenot_vs_gen_corr.pdf"
## pthr = 0.001 # threshold of nominal p-value based on FDR < 0.1


## 1. get genetic correlations -----------
# Starting here analysis
files = list.files(bi_est_dir, pattern = "_estNste.Rdata")

# 1. Parsing data ---------
stopifnot(length(files) > 0)
resVCs = lapply(files, function(f) {
  load(file.path(bi_est_dir,f))
  cat("doing file ", f, "\n")
  if(exists("res")){ 
    VCs = res$VCs
  }else{
    VCs = VCs
  }
  
  if(! "taskID" %in% colnames(VCs)){
    # TODO: this will be removed eventually
    # for now keep just because CFW do not have taskID col, because they where collected before modifying the script
    VCs[,"taskID"] = paste(VCs[,"trait1"], VCs[,"trait2"], sep='_') 
    order_cols = c("taskID",colnames(VCs))
    VCs = VCs[,order_cols]
  }
  if(any(duplicated( VCs[,"trait2"])) ) stop("pbm with duplicated trait2, move it to trait1 or comment rownames, might give a problem later")
  rownames(VCs) = VCs[,"trait2"] # this will raise a problem if trait 2 is duplicated; shouldn't be
  
  return(VCs)
})
names(resVCs) = gsub(paste0("_corr_.*"), "", gsub("_estNste.Rdata", "", files) ) #gsub("_estNste.Rdata", "", files)
#lapply(resVCs, head, 2)
#str(resVCs)

# Defining order of categories 
if (pop=="CFW mice"){
  dict_ord = c("immunology", "blood", "physiology", "anthropometric", "bones", "brain", "other", "behaviour") # behaviours last
}else if (pop == "HS mice"){
  dict_ord = c("immunology", "blood", "glucose", "steroid","physiology", "anthropometric", "brain", "other", "behaviour") # behaviours last
}else{
  stop("pb with population")
}

ormacro = read.csv(macro_file, header = T)
ormacro = ormacro[order(factor(ormacro$category, levels=dict_ord)), ]
category = setNames( ormacro$phenotype_ID, ormacro$category) 
macropheno = setNames( ormacro$phenotype_ID, ormacro$macrophenotype)

### specific CFW
if (pop =="CFW"){
  replacements = c("Adrenals" = "AdrenalWeight", 
                   "Bioch" = "Biochemistry", 
                   "WH" = "EarPunch",
                   "Haem" = "Haematology", 
                   "FC" = "Cue")
  names(macropheno) = sapply(names(macropheno), 
                             function(x, y) ifelse(x %in% names(y), return(y[x]), x), 
                             replacements, USE.NAMES = F)
}

# Add macropheno to resVCs
resRows = sapply(resVCs, function(x) nrow(x))
if (length( unique(resRows) ) > 1){
  # Align all to the min rows 
  refRows = rownames(resVCs [[which.min(resRows)]])
  resVCs = lapply(resVCs, function(x) x[refRows,])
}
sapply(resVCs, function(x) all(rownames(x) == rownames(resVCs[[1]])) )

# Adding corresponding macropheno 1 and macropheno 2
resVCs = lapply(resVCs, function(VCs){
  mt = match(VCs[,"trait1"], macropheno)
  VCs[,"macro1"] = names(macropheno)[mt]
  mt2 = match(VCs[,"trait2"], macropheno)
  VCs[,"macro2"] = names(macropheno)[mt2]
  
  # all df are going to have the same order
  #VCs = VCs[order(factor(VCs[,"macro2"], levels=names(dict))),]
  VCs = VCs[order(factor(VCs[,"macro2"], levels=unique(names(macropheno)))),]
  return(VCs)
})
# Checking all rownames are the same before cbind - i.e. they all have the same order
all(sapply(resVCs, function(x) all(rownames(x) == rownames(resVCs[[1]])) ))

# add category to resVCs
resRows = sapply(resVCs, function(x) nrow(x))
if (length( unique(resRows) ) > 1){
  # Align all to the min rows 
  refRows = rownames(resVCs [[which.min(resRows)]])
  resVCs = lapply(resVCs, function(x) x[refRows,])
}
sapply(resVCs, function(x) all(rownames(x) == rownames(resVCs[[1]])) )

# Adding corresponding category 1 and category 2
resVCs = lapply(resVCs, function(VCs){
  mt = match(VCs[,"trait1"], category)
  VCs[,"category1"] = names(category)[mt]
  mt2 = match(VCs[,"trait2"], category)
  VCs[,"category2"] = names(category)[mt2]
  
  # all df are going to have the same order
  #VCs = VCs[order(factor(VCs[,"macro2"], levels=names(dict))),]
  VCs = VCs[order(factor(VCs[,"category2"], levels=unique(names(category)))),]
  return(VCs)
})
# Checking all rownames are the same before cbind - i.e. they all have the same order
all(sapply(resVCs, function(x) all(rownames(x) == rownames(resVCs[[1]])) ))


# Ordering elements in list in function of order category of phenotypes
mt = match(names(resVCs), category)
mt2 = match(names(resVCs), macropheno)
resVCs = resVCs[order(factor(names(category)[mt], levels=unique(names(category))), 
                      factor(names(macropheno)[mt2], levels=unique(names(macropheno)))) ]

VCs.mat = do.call(cbind, resVCs) # matrix with all results in cbind
# Order rows according to category
VCs.mat = VCs.mat[category[match(rownames(VCs.mat), category)],] 


# getting names of trait1 and trait2 from VC matrix - going to need this to order and filter matrices
trait1 = gsub(".trait1", "", grep("trait1", colnames(VCs.mat), value = T))
trait2 = rownames(VCs.mat)
print("dim VCs.mat: ")
print(dim(VCs.mat))
trait1N2 = unique(c(trait1, trait2)[order(factor(c(trait1, trait2), levels=macropheno))])

corr_Ad2s1 = select_col(VCs.mat, "corr_Ad2s1")
colnames(corr_Ad2s1) = sub(".corr_Ad2s1", "", colnames(corr_Ad2s1))
corr_Ad2s1 = t(corr_Ad2s1)

corr_Ad1d2 = select_col(VCs.mat, "corr_Ad1d2")
colnames(corr_Ad1d2) = sub(".corr_Ad1d2", "", colnames(corr_Ad1d2))
corr_Ad1d2 = t(corr_Ad1d2)

## get phenotypic correlation
if(is.null(h5file)) stop("Need h5file to get phenotypic measurements")
if(is.null(phenoV)) stop("Need pheno version to get phenotypic measurements")
if(is.null(cageV)) stop("Need cage version to get cage assignment")

## 2. get phenos -----------
pheno = h5read(h5file, paste0("/phenotypes/",phenoV,"/matrix")); 
colnames(pheno) = h5read(h5file, paste0("/phenotypes/",phenoV,"/col_header/phenotype_ID"));
rownames(pheno) = h5read(h5file, paste0("/phenotypes/",phenoV,"/row_header/sample_ID"));
pheno[which(pheno == -999)] = NA
pheno = as.data.frame(pheno[!apply(pheno, 1, function(x) all(is.na(x))),])

cages = h5read(h5file, paste0("/cages/",cageV,"/array"))
names(cages) = h5read(h5file, paste0("/cages/",cageV,"/sample_ID"))

## 3. Phenotypic correlation focal trait - focal trait ---------
## Function getting the correlation between two phenotypes, either pearson or spearman, 
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

## 4. Phenotypic correlation with average of cage mates ---------
## Function to get phenotypic correlations matrix with average of cagemate on phenos2
cage_mate_cor_matrix = function(data, cages,
                                phenos1 = colnames(data),
                                phenos2 = colnames(data),
                                method = c("auto", "pearson", "spearman", "kendall"),
                                alpha_norm = 0.05, min_n_kendall = 20,
                                p.adjust.method = "fdr") {
  method = match.arg(method)
  
  ids     = intersect(rownames(data), names(cages))
  cage_of = cages[ids]
  X = as.matrix(data[ids, phenos1, drop = FALSE])   # focal phenotypes
  Y = as.matrix(data[ids, phenos2, drop = FALSE])   # phenotypes averaged over mates
  
  # cage mates (excluding self), computed once
  mates_list = lapply(seq_along(ids), function(i) {
    ids[cage_of == cage_of[i] & ids != ids[i]]
  })
  names(mates_list) = ids
  
  # matrix of cage-mate means: rows = individuals, cols = phenos2
  M = do.call(rbind, lapply(mates_list, function(m) {
    if (length(m) == 0) return(rep(NA_real_, length(phenos2)))
    colMeans(Y[m, , drop = FALSE], na.rm = TRUE)
  }))
  M[is.nan(M)] = NA_real_
  dimnames(M) = list(ids, phenos2)
  
  # correlation for one pair of vectors
  cor_pair = function(x, m) {
    ok = !is.na(x) & !is.na(m)
    xx = x[ok]; mm = m[ok]; n = length(xx)
    out = list(cor = NA_real_, p = NA_real_, n = n, method = NA_character_)
    if (n < 3 || sd(xx) == 0 || sd(mm) == 0) return(out)
    
    chosen = method
    if (method == "auto") {
      normal = if (n <= 5000) {
        all(c(shapiro.test(xx)$p.value, shapiro.test(mm)$p.value) > alpha_norm)
      } else TRUE
      if (normal) {
        chosen = "pearson"
      } else {
        tie_frac = 1 - min(length(unique(xx)), length(unique(mm))) / n
        chosen = if (n < min_n_kendall || tie_frac > 0.3) "kendall" else "spearman"
      }
    }
    res = tryCatch(cor.test(xx, mm, method = chosen, exact = FALSE),
                   error = function(e) NULL)
    if (is.null(res)) return(out)
    list(cor = unname(res$estimate), p = res$p.value, n = n, method = chosen)
  }
  
  dn = list(phenos1, phenos2)
  cor_mat = p_mat = n_mat = matrix(NA_real_, length(phenos1), length(phenos2), dimnames = dn)
  method_mat = matrix(NA_character_, length(phenos1), length(phenos2), dimnames = dn)
  
  for (i in seq_along(phenos1)) {
    for (j in seq_along(phenos2)) {
      r = cor_pair(X[, i], M[, j])
      cor_mat[i, j]    = r$cor
      p_mat[i, j]      = r$p
      n_mat[i, j]      = r$n
      method_mat[i, j] = r$method
    }
  }
  
  padj_mat = matrix(p.adjust(p_mat, method = p.adjust.method),
                    nrow = nrow(p_mat), dimnames = dn)
  
  # long format, convenient for filtering/sorting
  long = data.frame(
    pheno1 = rep(phenos1, times = length(phenos2)),
    pheno2 = rep(phenos2, each  = length(phenos1)),
    cor    = as.vector(cor_mat),
    p.value = as.vector(p_mat),
    p.adj  = as.vector(padj_mat),
    n      = as.vector(n_mat),
    method = as.vector(method_mat),
    stringsAsFactors = FALSE
  )
  
  list(cor = cor_mat, p.value = p_mat, p.adj = padj_mat,
       n = n_mat, method = method_mat, long = long,
       mate_means = M)
}
## get phenotypic matrix
# pheno IGE = pheno 1 = focal trait
# pheno DGE = pheno 2 = average cage mates trait
pheno_cm = cage_mate_cor_matrix(pheno, cages, phenos1 = rownames(corr_Ad2s1), phenos2 = colnames(corr_Ad2s1), method = "spearman")
pheno_mat_cm = pheno_cm$cor
dim(pheno_mat_cm)


## 5. plot -------
plot_cors = function(cormatx, cormaty, xlabi = "cormatx", ylabi = "cormaty", pointcolor = F, ...){
  cormatx = as.matrix(cormatx)
  cormaty = as.matrix(cormaty)
  plot(cormatx, cormaty, xlab = xlabi, ylab = ylabi, ylim = c(-1,1), xlim = c(-1,1), ...)
  abline(a = 0, b = 1, col = "grey", lty = 2)
  test = cor.test(cormatx, cormaty)
  plab = formatC(test$p.value, format = "e", digits = 2) 
  mtext(side = 1, line = 4, paste0("cor = ", round(test$estimate, 2), 
                         "; p-val = ", plab), cex = 0.8)
  if(pointcolor){
    points(cormatx[mask], cormaty[mask], 
           col = coolors[mask], pch = pch[mask])
  }
  return(test)
}
# Selecting p-val matrix to create mask for only significant ones, p-threshold is based on fdr significance
pmat = select_col(VCs.mat, "pv_chi2dof1")
colnames(pmat) = sub(".pv_chi2dof1", "", colnames(pmat))
pmat = t(pmat)
mask = pmat < pthr

dotcol = ifelse(pop == "HS mice", "darkorange", "#C51B7D") # as in phenot pipeline c("#E1A840","#8F407F") ; "#E7298AFF"
coolors = matrix("black", nrow = nrow(corr_Ad2s1), ncol= ncol(corr_Ad2s1), dimnames = dimnames(corr_Ad2s1))
coolors[mask] = adjustcolor(dotcol, alpha.f = 0.8)
pch = matrix(1, nrow = nrow(corr_Ad2s1), ncol= ncol(corr_Ad2s1), dimnames = dimnames(corr_Ad2s1))
pch[mask] = 19

pdf(outplot, w = 10, h = 5.5)
par(mfrow = c(1,2), las = 1, cex.lab = 1.2)
plot_cors(corr_Ad2s1, pheno_mat_cm, pointcolor = T,
          xlabi = "IGE-DGE correlation", ylabi = "phenotypic correlation (cm average)", main = pop)
plot_cors(corr_Ad1d2, pheno_mat, pointcolor = F,
          xlabi = "DGE-DGE correlation", ylabi = "phenotypic correlation", main = pop)

# only selected
plot_cors(corr_Ad2s1[mask], pheno_mat_cm[mask], main = pop,
          xlabi = "IGE-DGE correlation", ylabi = "phenotypic correlation (cm average)", 
          col = adjustcolor(dotcol, alpha.f = 0.8), pch = 19)

dev.off()


