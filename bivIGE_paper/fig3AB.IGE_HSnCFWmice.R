# Plot to compare IGE in behavioural vs non-behav
# 2 boxplots for each mouse dataset comparing behav and non-behav in terms of var_As

ymax = 0.32 # Same as Santostefano et al. 2025

#### INPUT FOR HS MICE - comment input for CFW mice #####
pop="HSmice"; popmain="HS mice"
unires = "~/PRJs/HSmice/output/VDreal_2311/univariate/data_bcNcovariates_500/Andres_kinship_DGE_IGE_cageEffect_estNste.Rdata"
macro="~/PRJs/HSmice/output/dataset/macropheno_HSmice.csv"
colbycat = F; outfile = "./plot/HSmice/Fig3.IGE_bplot_behav_nonbeh.pdf"
# To colour by phenotypic category instead of by beh / non-beh uncomment below
#colbycat = T; outfile = "./plot/HSmice/IGE_bplot_behav_nonbeh_colbycat.pdf"


#### INPUT FOR CFW MICE - comment input for HS mice #####
# pop = "CFW"; popmain="CFW mice"
# unires = "~/PRJs/CFW/output/VDreal_2310/univariate/noBatch_500/pruned_dosages_include_DGE_IGE_cageEffect_estNste.Rdata"
# macro="~/PRJs/CFW/output/dataset/macropheno_CFW.csv"
# colbycat = F; outfile = "./plot/CFW/Fig3.IGE_bplot_behav_nonbeh.pdf"
# # To colour by phenotypic category instead of by beh / non-beh uncomment below
# #colbycat = T; outfile = "./plot/CFW/IGE_bplot_behav_nonbeh_colbycat.pdf"


#### defining order of categories ####
if (pop=="CFW"){
  dict_ord = c("immunology", "blood", "physiology", "anthropometric", "bones", "brain", "other", "behaviour") # behaviours last
}else if (pop == "HSmice"){
  dict_ord = c("immunology", "blood", "glucose", "steroid","physiology", "anthropometric", "brain", "other", "behaviour") # behaviours last
}else{
  stop("pb with population")
}

# Read file with info on macrophenotypes
ormacro = read.csv(macro, header = T)
ormacro = ormacro[order(factor(ormacro$category, levels=dict_ord)), ]
# Define categories as behaviours / non-behaviours
new_cat = c("immunology" = "non-behaviour",
            "blood" = "non-behaviour",
            "glucose" = "non-behaviour",
            "steroid" = "non-behaviour",
            "physiology" = "non-behaviour",
            "bones" = "non-behaviour",
            "brain" = "non-behaviour",
            "other" = "non-behaviour",
            "anthropometric" = "non-behaviour",
            "behaviour" = "non-social\nbehaviour"
)
# Add new category to macrophenotype df
ormacro[,"beh_nonb"] = new_cat[ormacro[,"category"]]
table(ormacro[,"beh_nonb"])
macropheno = setNames( ormacro$phenotype_ID, ormacro$macrophenotype)
cat = setNames( ormacro$phenotype_ID, ormacro$category)
newcat = setNames( ormacro$phenotype_ID, ormacro$`beh_nonb`)

### specific CFW - to align to HSmice
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

# Load results from univariate 
load(unires) 
# !!NB: loading object 'uni_VCs' as saved like that initially, if running latest version of the VD code, need to uncomment line below
#uni_VCs = res$VCs
uni_VCs = uni_VCs[macropheno,]
all(rownames(uni_VCs) == macropheno)
#uni_VCs[,"macropheno"] = names(macropheno)
uni_VCs[,"category"] = names(cat)
uni_VCs[,"beh_nonb"] = names(newcat)

# Define x positions for points
at.dict = sapply(uni_VCs[,"beh_nonb"], function(x) ifelse(x =="non-social\nbehaviour", 2, 1))

if(colbycat){
  # colour all categories
  colcat = c("immunology"            = "#53585F",
             "blood"                 = "#5D6D84",
             "glucose"               = "#8DA0CB",
             "steroid"               = "#9EC2C9",
             "physiology"            = "#8ABFAA",
             "anthropometric"        = "#A4CBA6",
             "bones"                 = "#88A77D",
             "brain"                 = "#5C8A72",
             "other"                 = "#004242",
             "non-social\nbehaviour" = "#B894B1")
  coolors = colcat[uni_VCs[,"category"]]
  
}else{
  # two colours
  coolors = sapply(at.dict, function(x) ifelse(x == 2, "#B894B1", "#70928D"))
}


######## Plotting IGE ############
pdf(outfile, h= 6, w=4)
par(mar = c(6.1, 5.1, 3.1, 1.1))
farmula = as.formula(prop_As1 ~ beh_nonb) 
#boxplot(farmula) # simple boxplot
#ylimi = c(0, max(uni_VCs[,"prop_As1"])+0.025) # set limit in function of the dataset
ylimi = c(0 - 0.005, ymax + 0.005) # set limit in function of ymax defined at the beginning (if want to have same y-axis across datasets)
bp = boxplot(farmula, data = uni_VCs, 
             col = adjustcolor("white", 0), border=adjustcolor("white", 0), 
             ylim = ylimi,
             outline=T, ylab = "IGE", xlab = '', xaxt="n",
             main = popmain, drawRect=F, varwidth = TRUE, las = 1, cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.5)
# add ticks and ticks labels to x-axis
axis(1, at=unique(at.dict), labels=F, cex.axis=1.25)
axis(1, at=unique(at.dict), labels=unique(names(at.dict)), line = 0.5, cex.axis=1.2, tick=F, lwd=0)

# add points
set.seed(3)
points(x = jitter(unname(at.dict), factor=0.5), 
       y = uni_VCs[,"prop_As1"], 
       pch=16, cex=1, 
       col = coolors)
# add boxes
bp = boxplot(farmula, data = uni_VCs, 
             col = adjustcolor("white", alpha = 0.7),
             outline=F, 
             ylab = "", bty = "n",
             xlab = '', xaxt="n",yaxt="n",
             drawRect=F,
             varwidth = F,boxwex = 0.5,
             las = 1, cex.axis = 1.2, cex.main = 1.5,
             add=T)  

# using wilcox test - given our estimates are not normally distributed
# testing for non behaviour to be smaller than behaviour
x = uni_VCs[uni_VCs$beh_nonb == "non-social\nbehaviour", "prop_As1"]
y = uni_VCs[uni_VCs$beh_nonb == "non-behaviour", "prop_As1"]
wtest = wilcox.test(x, y, alternative = "greater") # testing for x, behav, being greater than y, non-behav --> non-sign

# testing for difference between beh and non-beh 
#wtest = wilcox.test(farmula, data = uni_VCs, exact = FALSE)

alpha = 0.05
t = ifelse(wtest$p.value < alpha, "*", "ns")

segments(1, max(uni_VCs[,"prop_As1"]) + 0.01, 2, max(uni_VCs[,"prop_As1"]) + 0.01)
text(1.5, max(uni_VCs[,"prop_As1"]) + 0.02, t, font=2, cex = 1.2)
dev.off()
