suppressMessages(library("rhdf5"))
suppressMessages(library("dendextend")) # used for clustering
suppressMessages(library("here"))

sourcefun= "./bivIGE_paper/Rfun/"

source(here(sourcefun, "select_col_VCsmat.R")) # `select_col`: function to select columns of a specific Variance component, for all phenotypes 

# functions needed for my heatmap
source(here(sourcefun, "corrplot_size.R")) # `corrplot_size`: function to get heatmap - based on corrplot and with option "size_vector" 
source(here(sourcefun, "heatmapCore.R")) # `corrplot_legend`: function to get heatmap with size based on pvalue - using option "size_vector" - takes corrplot_size
source(here(sourcefun, "cluster_for_heatmap.R")) # `` : functions to get clustering - to plot clustered heatmap 


# To plot results for CFW mice
## opt = list(pop="CFW", # population, I have CFW or HSmice
##            # directory with bivariate estimates 
##            biVCdir ="~/nf_PRJs/nf-CoreQuantGen/realdata/output/CFW/VD/bivariate/noBatch/pruned_dosages_include_DGE_IGE_cageEffect/",
##            # file with macrophenotypes - needed with biVCdir
##            macro="~/PRJs/CFW/output/dataset/macropheno_CFW.csv", 
##            # dir where to save plot # can do default = "./"
##            out="./plot/CFW/") 

# To plot results for HS mice
opt = list(pop = "HSmice", # population
           # directory with bivariate estimates 
           biVCdir = "~/nf_PRJs/nf-CoreQuantGen/realdata/output/HSmice/VD/bivariate/data_bcNcovariates/Andres_kinship_None_DGE_IGE_cageEffect/",
           # file with macrophenotypes - needed with biVCdir
           macro = "~/PRJs/HSmice/output/dataset/macropheno_HSmice.csv", 
           # dir where to save plot # can do default = "./"
           out = "./plot/HSmice/") 


# 0. Storing options ---------
pop = opt$pop 
bi_est_dir = opt$biVCdir
macro_file = opt$macro
coroi = "corr_Ad2s1"
pval_type = "LRT"
adj_toplot = "fdr"
psig = "0.1"

files = list.files(bi_est_dir, pattern = "_estNste.Rdata")

# 1. Parsing data ---------
stopifnot(length(files) > 0)
resVCs <- lapply(files, function(f) {
  load(file.path(bi_est_dir,f))
  cat("doing file ", f, "\n")
  if(exists("res")){ 
    VCs = res$VCs
    #VCs0 = res$VCs0
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
if (pop=="CFW"){
  dict_ord = c("immunology", "blood", "physiology", "anthropometric", "bones", "brain", "other", "behaviour") # behaviours last
}else if (pop == "HSmice"){
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
  replacements <- c("Adrenals" = "AdrenalWeight", 
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
resVCs <- lapply(resVCs, function(VCs){
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
resVCs <- lapply(resVCs, function(VCs){
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


# Getting names of trait1 and trait2 from VC matrix - going to need this to order and filter matrices
trait1 = gsub(".trait1", "", grep("trait1", colnames(VCs.mat), value = T))
trait2 = rownames(VCs.mat)
print("dim VCs.mat: ")
print(dim(VCs.mat))
trait1N2 = unique(c(trait1, trait2)[order(factor(c(trait1, trait2), levels=macropheno))])

# 2. Getting cor_mat and p_mat - calculate pvalue if necessary 
#    Levels of significance to star in the plot
siglev = as.numeric(ifelse(is.na(psig), 0.1, psig)) # this is for the * to display


# 3. cor_mat and p_mat FOR ESTIMATED CORS -----------
#  Get ste name depending on corr name - and name for pvalue
ste_oi = gsub("corr", "STE", coroi) #"STE_Ad2s1", #"STE_Ad1d2" # "STE_As1s2" #
p_oi = gsub("STE_", "pnom_", ste_oi)
p_oi.fdr = gsub("STE_", "fdr_", ste_oi) # could add the bonferroni, not for now

#    Calculating pvalues
#   This is for now because have the null model only for corr_Ad2s1
#stopifnot(coroi=="corr_Ad2s1")

#    Extracting pvalue from model estimate
p_nom = as.matrix(select_col(VCs.mat, "pv_chi2dof1")) # 
colnames(p_nom) = gsub("pv_chi2dof1",coroi,colnames(p_nom))

## the whole matrix together
cat("Calculating FDR on the full matrix\n")
pval_flat = c(as.matrix(p_nom))
## Adjust as a vector 
p_adj.2 <- p.adjust(pval_flat, method = "fdr")
## Reshape the adjusted p-values back into a matrix
p_adj.2mx <- matrix(p_adj.2, nrow = nrow(p_nom), dimnames = dimnames(p_nom))

p_fdr = p_adj.2mx

stopifnot(all.equal(rownames(p_nom), rownames(p_fdr)), 
          all.equal(colnames(p_nom), colnames(p_fdr)))


#  Extracting the matrix for the corr of interest
cor_mat <- select_col(VCs.mat, coroi) 
se_mat <- select_col(VCs.mat, ste_oi) 
stopifnot(all(rownames(cor_mat) == rownames(se_mat)), 
          all(gsub("corr_","",colnames(cor_mat)) == gsub("STE_", "", colnames(se_mat))))
# plot(unlist(as.vector(cor_mat)), unlist(as.vector(se_mat))) # checking that there is no correlation between corr value and ste value

p_mat = p_fdr

# Ordering so that order the same as phenotypic results
mot = na.omit(match(trait1N2, rownames(cor_mat)))
cor_mat = cor_mat[mot,]
se_mat = se_mat[mot,]
p_mat = p_mat[mot,]

# Checking the distribution of estimates
#hist(unlist(c(cor_mat)), main=paste0("distribution of estimates - ", pop), xlab=coroi)
#box(bty="l")

# At this point I have
#   dim(VCs.mat)
#   dim(cor_mat)
#   dim(p_nom)
#   dim(p_mat)
#   siglev
#   adj_toplot

#    Define the title - based on the siglevels and padjust
ttl = vector(length=0L)
for (i in seq_along(siglev)){
  ttl = c(ttl, paste0("p(",adj_toplot,")<", rev(siglev)[i],  ": ", paste(rep("*", i), collapse="")))
}
plot.title = paste0(coroi, " - significance: ", paste(ttl, collapse = " , "), " - ", pop)
cat(plot.title, "\n")

# Making sure elements in cor_mat map to element in p_mat # 130 x 35 or 100 x 42 - rows DGE phenos, cols IGE phenos
stopifnot(all(rownames(cor_mat)==rownames(p_mat) ),  
          all(colnames(cor_mat)==colnames(p_mat) ) )
# NB: to map in corplot they need to have the same colnames

mat.toplot = t(cor_mat[trait2, paste0(trait1, ".", coroi)])
p_mat.toplot = t(p_mat[trait2, paste0(trait1, ".", coroi)])


## clustering rows and columns - euclidean
tree_row = cluster.mat(p_mat.toplot, distance = "euclidean", method = "complete")
tree_col = cluster.mat(t(p_mat.toplot), distance = "euclidean", method = "complete")
## extracting the order from the clustering 
row_ord = tree_row$order
col_ord = tree_col$order


## Reordering for HS mice
if(pop=="HSmice"){
  # Convert row tree to dendrogram
  dend_row <- as.dendrogram(tree_row)
  #par(mar=c(10.1,4.1,4.1,2.1))
  #nodePar = list(lab.cex = 0.5, pch=c(NA,NA))
  #plot(dend_row, nodePar = nodePar)
  # Reorder rows 
  dend_row[[2]][[2]] <-  rev(dend_row[[2]][[2]])
  dend_row[[2]][[2]][[1]] <-  rev(dend_row[[2]][[2]][[1]])
  dend_row[[2]][[2]][[1]][[2]][[1]] <-  rev(dend_row[[2]][[2]][[1]][[2]][[1]])
  
  row_ord = order.dendrogram(dend_row)
  # saving new row tree
  tree_row = as.hclust(dend_row)
  
  # Convert col tree to dendrogram
  dend_col <- as.dendrogram(tree_col)
  #plot(dend_col, nodePar = nodePar)
  # Reorder cols 
  dend_col[[2]][[2]] <-  rev(dend_col[[2]][[2]])
  dend_col[[2]][[2]][[1]][[1]][[1]][[2]] <- rev(dend_col[[2]][[2]][[1]][[1]][[1]][[2]])
  dend_col[[2]][[2]][[1]][[2]] <- rev(dend_col[[2]][[2]][[1]][[2]])
  dend_col[[2]][[1]] <- rev(dend_col[[2]][[1]])
  dend_col[[2]][[1]][[2]] <- rev(dend_col[[2]][[1]][[2]])
  
  col_ord = order.dendrogram(dend_col)
  # saving new col tree
  tree_col = as.hclust(dend_col)
  
}else if(pop == "CFW"){
  # Convert row tree to dendrogram
  dend_row <- as.dendrogram(tree_row)
  # Reorder rows
  #dend_row[[1]] <- rev(dend_row[[1]])
  dend_row[[1]][[1]][[2]] <- rev(dend_row[[1]][[1]][[2]])
  dend_row[[1]][[1]][[2]][[1]] <- rev(dend_row[[1]][[1]][[2]][[1]])
  dend_row[[1]][[1]][[2]][[1]][[2]] <- rev(dend_row[[1]][[1]][[2]][[1]][[2]] )
  dend_row[[1]][[2]] <- rev(dend_row[[1]][[2]])
  dend_row[[2]] <- rev(dend_row[[2]])
  dend_row[[2]][[1]] <- rev(dend_row[[2]][[1]])
  dend_row[[2]][[1]][[2]][[1]] <- rev(dend_row[[2]][[1]][[2]][[1]])
  dend_row[[2]][[1]][[2]][[1]][[1]] <- rev(dend_row[[2]][[1]][[2]][[1]][[1]])
  dend_row[[2]][[2]] <- rev(dend_row[[2]][[2]])
  
  row_ord = order.dendrogram(dend_row)
  # saving new row tree
  tree_row = as.hclust(dend_row)
  
  # Convert col tree to dendrogram
  dend_col <- as.dendrogram(tree_col)
  # Reorder
  dend_col <- rev(dend_col)
  dend_col[[1]] <- rev(dend_col[[1]])
  dend_col[[1]][[2]][[2]] <- rev(dend_col[[1]][[2]][[2]])
  dend_col[[1]][[2]][[2]][[1]] <- rev(dend_col[[1]][[2]][[2]][[1]])
  dend_col[[2]][[2]] <- rev(dend_col[[2]][[2]])
  dend_col[[2]][[2]][[1]][[2]] <- rev(dend_col[[2]][[2]][[1]][[2]])
  dend_col[[2]][[2]][[1]][[2]][[1]] <- rev(dend_col[[2]][[2]][[1]][[2]][[1]])
  dend_col[[2]][[2]][[1]][[2]][[1]][[2]] <- rev(dend_col[[2]][[2]][[1]][[2]][[1]][[2]])
  dend_col[[2]][[2]][[2]][[1]][[1]] <- rev(dend_col[[2]][[2]][[2]][[1]][[1]])
  dend_col[[2]][[2]][[2]][[1]][[1]][[1]] <- rev(dend_col[[2]][[2]][[2]][[1]][[1]][[1]])
  dend_col[[2]][[2]][[2]][[1]][[1]][[1]][[2]][[2]] <- rev(dend_col[[2]][[2]][[2]][[1]][[1]][[1]][[2]][[2]])
  
  col_ord = order.dendrogram(dend_col)
  # saving new col tree
  tree_col = as.hclust(dend_col)
  
}

mat.toplot = mat.toplot[row_ord, col_ord]
p_mat.toplot = p_mat.toplot[row_ord,col_ord]

# 6. Create OUTPUT DIR depending on SET ----
outDir = file.path(opt$out)
dir.create(outDir, showWarnings = F)
cat("created output directory:", outDir, "\n")

# C. Heatmap - size as pvalues --------
# select the h - height, w - width, and tl.cex - label size
h = 12
w = 30
if(pop=="HSmice"){
  tl.cex = 1.4
}else if(pop=="CFW"){
  tl.cex = 1.2
}
## Text size when plot with legend

# This is to remove the coroi from the labels (which I don't like because easy to lose track but looks better in the plot)
rownames(mat.toplot) = rownames(p_mat.toplot) = gsub(paste0(".",coroi), "", rownames(mat.toplot)) # 

outpdf = file.path(outDir, paste0("fig5.",coroi,"_heatmap.pdf")); cat("saving plot in ", outpdf, "\n")

# Heatmap - with clustering based on corr value --------
pdf(outpdf, h=h, w =w)
corrplot_legend(mat.toplot, p_mat.toplot,
                siglev = siglev, plot.title = plot.title, tl.cex = tl.cex, cl.cex = 1, 
                pch.cex = 1.7, pch.col = "yellow")
dev.off()

library("ggsci")
# Plot tree of the clustering
#mycolors = pal_npg("nrc")(10) #c("black", "#DF536B", "#61D04F", "#2297E6", "#28E2E5", "#CD0BBC", "#F5C710", "gray62", "darkviolet", "darkorange3") # = 1:8 ; + other two
mycolors = c("#e40203","#ff8b00", "#feed01", "#007f24", "#004dff", "#760789")

outpdf = file.path(outDir, paste0("fig5.",coroi,"_tree.pdf")); cat("saving plot in ", outpdf, "\n")
h.t=6; w.t = 15
pdf(outpdf, h=h.t, w =w.t)
k = 6
plot.tree(tree_col, k = k, subtitle= "tree cols - DGE2", mycolors, cex = 0.5)
plot.tree(tree_row, k = k, subtitle= "tree rows - IGE1", rev(mycolors), cex = 0.8)

dev.off()
