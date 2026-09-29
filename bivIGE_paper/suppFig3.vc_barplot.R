#############################################################
############  Script to plot Supp Fig. 3   ##################
# Plot to compare magnitude DGE and IGE in the two datasets #
#############################################################

maxbr=0.7

prop = "absolute"; ylabi = "number of traits"
maxy=135 # absolute n - DGE as IGE

var="prop_Ad1"; xlabi="DGE"
var="prop_As1"; xlabi="IGE"

## To plot HSmice - comment lines CFW
#pop = "HSmice"; popmain="HS mice"
#unires = "./data/HSmice/VD/univariate/data_bcNcovariates_Andres_kinship_DGE_IGE_cageEffect_estNste.Rdata"
#macro = "./data/HSmice/dataset/macropheno_HSmice.csv"
#outfile = paste0("./plot/HSmice/SFig3.VC_barplot_",xlabi,"_",prop,".pdf")

## To plot CFW - comment lines HSmice
pop = "CFW"; popmain="CFW mice"
unires = "./data/CFW/VD/univariate/noBatch_pruned_dosages_include_DGE_IGE_cageEffect_estNste.Rdata"
macro = "./data/CFW/dataset/macropheno_CFW.csv"
outfile = paste0("./plot/CFW/SFig3.VC_barplot_",xlabi,"_",prop,".pdf")



#### defining order of categories ####
if (pop=="CFW"){
  dict_ord = c("immunology", "blood", "physiology", "anthropometric", "bones", "brain", "other", "behaviour") # behaviours last
}else if (pop == "HSmice"){
  dict_ord = c("immunology", "blood", "glucose", "steroid","physiology", "anthropometric", "brain", "other", "behaviour") # behaviours last
}else{
  stop("pb with population")
}

ormacro = read.csv(macro, header = T)
ormacro = ormacro[order(factor(ormacro$category, levels=dict_ord)), ]
table(ormacro[,"category"])
macropheno = setNames( ormacro$phenotype_ID, ormacro$macrophenotype)
cat = setNames( ormacro$phenotype_ID, ormacro$category)

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

load(unires)# loads uni_VCs
#uni_VCs_safe = uni_VCs
uni_VCs = uni_VCs[macropheno,]
all(rownames(uni_VCs) == macropheno)
#uni_VCs[,"macropheno"] = names(macropheno)
uni_VCs[,"category"] = names(cat)

#### mean and prevalence of IGE
## mean(uni_VCs$prop_As1)
## 
## sum(uni_VCs[which(uni_VCs$category == "behaviour"),"prop_As1"] > 0.05) # behav pheno IGE > 0.05
## nrow(uni_VCs[which(uni_VCs$category == "behaviour"),]) # tot behav pheno
## mean(uni_VCs[which(uni_VCs$category == "behaviour"),"prop_As1"]) # mean behav pheno
## 
## sum(uni_VCs[which(uni_VCs$category != "behaviour"),"prop_As1"] > 0.05) # non-behav pheno IGE > 0.05
## nrow(uni_VCs[which(uni_VCs$category != "behaviour"),]) # tot non-behav pheno
## mean(uni_VCs[which(uni_VCs$category != "behaviour"),"prop_As1"]) # mean non-behav pheno

coolors = c("immunology"     = "#53585F", #"#6E8790",
            "blood"          = "#5D6D84", #"#7F8F9D",
            "glucose"        = "#8DA0CB",#"#83989F", "#729EAB",#"#7E9B45", 
            "steroid"        = "#9EC2C9",
            "physiology"     = "#8ABFAA",#"#5BAF95", #"#66C2A5",
            "anthropometric" = "#A4CBA6",
            "bones"          = "#88A77D",
            "brain"          = "#5C8A72", #"#8ABFAA", "#539D88","#5BAF95",
            "other"          = "#004242",#"#55A868",
            "behaviour"      = "#B894B1") #"#A05D79")


##### Plotting vertical barplot ##########
# Define function to obtain counts for each set of phenotypes to plot
#     subdf: sub-dataframe with set of phenotypes of interest
#     var: proportional variance to plot - ideally could do this for other 'prop_'
#     br: breaks in histogram

prepare_hist = function(subdf, var = "prop_Ad1", br, prop = "relative"){
  hist_values = hist(subdf[, var], breaks = br, plot = F)
  if(prop=="relative"){ #relative to subdf
    counts = hist_values$counts/sum(nrow(subdf))
  }else if(is.numeric(prop)){
    counts = hist_values$counts/prop
  }else if(prop=="absolute"){
    counts = hist_values$counts
  }
  return(counts)
}

# defining breaks; if want to compare different subset in the same plot, need breaks to be all the same
br = seq(0, maxbr, by = 0.05) # flexible, defined at the beginning

# creating dataframe to plot
toplot = t(sapply(names(coolors), function(c){prepare_hist(uni_VCs[uni_VCs[,"category"]==c, ], 
                                                           var, 
                                                           br, 
                                                           prop=prop)} ))
colnames(toplot) = br[-1]
# removing lines with all NA values
toplot = toplot[!apply(is.na(toplot), 1, all), ]

par(mar = c(5.1, 5.1, 2.1, 2.1)) # default: c(5.1, 4.1, 4.1, 2.1)

# Set ylim 
ylimi=c(0,maxy) # flexible, defined at the beginning

cat("saving file to ", outfile, "\n")
pdf(outfile, w= 10, h = 7)

# plot
bp <- barplot(toplot,
        beside = F, names.arg = rep("", ncol(toplot)), #colnames(toplot), 
        col = coolors[rownames(toplot)], las = 1, 
        main=popmain, cex.main=2,
        ylim = ylimi,  
        horiz=F,
        cex.names = 1.25, cex.axis = 1.25, cex.lab=1.4,
        ylab = ylabi, 
        xlab=xlabi)

# Adding axis with values in the middle of the columns
xpos = bp + c(diff(bp)/2, (diff(bp)/2)[length(bp)-1])
axis(1, at = c(0, xpos), 
     labels = c("0.0",colnames(toplot)),
     cex.axis = 1.25, tick = F, lwd=0)

# subtitle with tot number of traits
mtext(side=3, line=0, text= paste0("(tot traits = ", nrow(uni_VCs), ")"), cex=1.4)

label_positions <- bp[-length(bp)] + diff(bp)/2

# Add a semi-transparent gray bar on top of the first bar
# Get the height of the first bar
first_bar_height <- sum(toplot[, 1])  # Sum of all values in the first stacked bar

# Add a semi-transparent gray rectangle over the first bar - VC <0.05
rect(bp[1] - 0.5,         # Left side of the bar
     0,                   # Bottom of the bar
     bp[1] + 0.5,         # Right side of the bar
     first_bar_height,    # Height of the bar
     col = rgb(1, 1, 1, 0.6),  # Semi-transparent gray
     border = rgb(1, 1, 1, 0.6))  

# legend
legend (x = 'topright', 
        legend = rownames(toplot), 
        fill = coolors[rownames(toplot)], border = NA, cex =1.2, bty="o", inset=c(0.04,0)) 
# inset so that box of legend aligns to last bin - with margins par(mar = c(5.1, 5.1, 2.1, 2.1))  and width = 8
dev.off()


