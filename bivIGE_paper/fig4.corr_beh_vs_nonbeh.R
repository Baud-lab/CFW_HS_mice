### HSmice - comment lines for CFW ####
# pop = "HSmice"; pop_main = "HS mice"
# bi_est_dir = "~/nf_PRJs/nf-CoreQuantGen/realdata/output/HSmice/VD/bivariate/data_bcNcovariates/Andres_kinship_None_DGE_IGE_cageEffect/"
# macro_file = "~/PRJs/HSmice/output/dataset/macropheno_HSmice.csv"
# outfile = "./plot/HSmice/Fig4.corr_beh_vs_nonb_merge.pdf"

### CFW - comment lines for HSmice ####
pop = "CFW"; pop_main = "CFW mice"
bi_est_dir = "~/nf_PRJs/nf-CoreQuantGen/realdata/output/CFW/VD/bivariate/noBatch/pruned_dosages_include_DGE_IGE_cageEffect/"
macro_file = "~/PRJs/CFW/output/dataset/macropheno_CFW.csv"
outfile = "./plot/CFW/Fig4.corr_beh_vs_nonb_merge.pdf"


# Listing files in output directory
files = list.files(bi_est_dir, pattern="_estNste.Rdata", full.names = F) #[1:9]

# Getting results - output a list of VCs
#f=files[1]
stopifnot(length(files) > 0)
resVCs <- lapply(files, function(f) {
  load(file.path(bi_est_dir,f))
  VCs = res$VCs
  #VCs0 = res$VCs0
  
  if(any(duplicated( VCs[,"trait2"])) ) stop("pbm with duplicated trait2, move it to trait1 or comment rownames, might give a problem later")
  rownames(VCs) = VCs[,"trait2"] # this will raise a problem if trait 2 is duplicated; shouldn't be
  
  return(VCs)
})
# Set names of list as basename of files
names(resVCs) = gsub(paste0("_corr_.*"), "", gsub("_estNste.Rdata", "", files) )
resVCs = Filter(Negate(is.null), resVCs)
# See results 
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
# Change categories
new_cat = c("immunology" = "non-behaviour",
            "blood" = "non-behaviour",
            "glucose" = "non-behaviour",
            "steroid" = "non-behaviour",
            "physiology" = "non-behaviour",
            "bones" = "non-behaviour",
            "brain" = "non-behaviour",
            "other" = "non-behaviour",
            "anthropometric" = "non-behaviour",
            "behaviour" = "behaviour"
)
# Add the categories to macrophenotype df
ormacro[,"category"] = new_cat[ormacro[,"category"]]
category = setNames( ormacro$phenotype_ID, ormacro$category) 
macropheno = setNames( ormacro$phenotype_ID, ormacro$macrophenotype)

# specific CFW - set names as in HSmice
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
  VCs = VCs[order(factor(VCs[,"macro2"], levels=unique(names(macropheno)))),]
  return(VCs)
})
# Checking all rownames are the same before cbind - i.e. they all have the same order
all(sapply(resVCs, function(x) all(rownames(x) == rownames(resVCs[[1]])) ))


# Add category to resVCs
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

# Slim results
colkeep = c("trait1", "trait2", "category1", "category2", "prop_Ad1", "prop_Ad2", "prop_As1", "prop_As2", "corr_Ad1s1", "corr_Ad1s2", "corr_Ad2s1", "STE_Ad1s1", "STE_Ad2s1", "STE_Ad1s2", "pv_chi2dof1")
slim_res = lapply(resVCs, function(x){x[,colkeep]})

# See dataframe
lapply(slim_res, head, 2)


# Select the same trait 
# plot on one line
# select the non-same trait
# plot on the other line

# Gather all results in the same df
res = do.call(rbind,  slim_res)
dim(res)

set.seed(11)
# Define positions on y-axis (with jitter)
res[, "y_pos"] = jitter(sapply(res$category1, function(x) ifelse(x == "behaviour", 3, 1.5)), amount = 0.05) # DGE-behaviour
res[res$category2 == "non-behaviour", "y_pos"] = res[res$category2 == "non-behaviour", "y_pos"] - 0.5 # DGE-non-behaviour

# Get mean and sd to plot big dot and lines
mean_sd = vapply(unique(res$category1), function(x){
  sel = res[res[,"category1"] == x,]
  beh = sel["category2"] == "behaviour"
  mean_beh = mean(abs(sel[beh, "corr_Ad2s1"]))
  sd_beh = sd(abs(sel[beh, "corr_Ad2s1"]))
  
  mean_nonb = mean(abs(sel[!beh, "corr_Ad2s1"]))
  sd_nonb = sd(abs(sel[!beh, "corr_Ad2s1"]))
  
  # Check distribution of estimates
  #   hist(abs(sel[beh, "corr_Ad2s1"]))
  #   hist(abs(sel[!beh, "corr_Ad2s1"]))
  # Testing if correlations from behaviour(DGE) are stronger than from non-behaviour(DGE)
  pval = wilcox.test(abs(sel[beh, "corr_Ad2s1"]), abs(sel[!beh, "corr_Ad2s1"]), alternative = "greater")$p.value
  return(c(mean_beh = mean_beh, sd_beh = sd_beh, mean_nonb = mean_nonb, sd_nonb=sd_nonb, pval=pval))
}, FUN.VALUE = c(mean_beh=0, sd_beh=0, mean_nonb=0, sd_nonb=0, pval=NA))

# Show p-values
print(mean_sd["pval",])

# Define position of ticks and labels on y-axis
y_pos = c("behaviour" = c("behaviour" = 3, "non-behaviour" = 2.5), 
          "non-behaviour" = c("behaviour" = 1.5, "non-behaviour" = 1))  # Reverse order for variables

coolors = c("behaviour"= "#B894B1", 
            "non-behaviour" = "#70928D")
#sapply(res$category2, function(x) ifelse(x == "behaviour", "#B894B1", "#70928D"))

# Define function to plot big dot and lines (mean±sd)
dot_lines = function(m, s, yp, col){
  lines(c(m-s, m+s), 
        c(yp, yp), 
        col = "black", 
        #col=col,
        lwd = 3)
  points(m, yp, 
         pch=21, col="black", cex=2.5,
         bg=col)
}

# Plot
pdf(outfile, h=5, w = 8)
par(mar=c(5.1, 15.1, 4.1, 1.1), cex.lab = 1.4, cex.axis=1.2)

plot(abs(res[,"corr_Ad2s1"]), res[,"y_pos"], pch = 16, cex = 1.2, col = adjustcolor("grey", alpha.f = 0.5), 
     yaxt="n",ylab="", xlab="IGE-DGE correlation", main = pop_main, cex.main = 1.5)
grid(ny=NA, nx=NULL)
axis(2, at=y_pos, labels = unlist(sapply(strsplit(names(y_pos), "[.]"), "[[", 2)), las = 1)
axis(2, at=c(2.75, 1.25), labels = unique(unlist(sapply(strsplit(names(y_pos), "[.]"), "[[", 1))), las = 1, 
     line = 7.5, tick=F, lwd=0)
axis(2, at=3.3, labels = "IGE", las = 1, font = 4, 
     line = 9, tick=F, lwd=0, xpd=T)
axis(2, at=3.3, labels = "DGE", las = 1, font = 4, 
     line = 1.5, tick=F, lwd=0, xpd=T)
dot_lines(mean_sd["mean_beh","behaviour"], 
          mean_sd["sd_beh","behaviour"], 
          y_pos["behaviour.behaviour"], 
          coolors["behaviour"])
dot_lines(mean_sd["mean_nonb","behaviour"], 
          mean_sd["sd_nonb","behaviour"], 
          y_pos["behaviour.non-behaviour"], 
          coolors["non-behaviour"])
dot_lines(mean_sd["mean_beh","non-behaviour"], 
          mean_sd["sd_beh","non-behaviour"], 
          y_pos["non-behaviour.behaviour"], 
          coolors["behaviour"])
dot_lines(mean_sd["mean_nonb","non-behaviour"], 
          mean_sd["sd_nonb","non-behaviour"], 
          y_pos["non-behaviour.non-behaviour"], 
          coolors["non-behaviour"])
dev.off()
