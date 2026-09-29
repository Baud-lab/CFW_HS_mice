#######################################################################
##################  Script to plot Supp Fig. 5   ######################
# Plot estimates from scrambling cage mates assignment (permutations) #
#######################################################################

## input output
fullfileHS = "./input/HSmice/VD/bivariate/data_bcNcovariates_Andres_kinship_None_DGE_IGE_cageEffect_corr_Ad1s2_alt_all.txt"
permfileHS = "./input/HSmice/permutations/bivariate/Andres_kinship_permALL_DGE_IGE_IEE_cageEffect_all_estNste.Rdata"

fullfileCFW = "./input/CFW/VD/bivariate/noBatch_pruned_dosages_include_DGE_IGE_cageEffect_corrAd1s2_alt_all.txt"
permfileCFW = "./input/CFW/permutations/bivariate/pruned_dosages_permALL_DGE_IGE_IEE_cageEffect_all_estNste.Rdata"

outpdf = "./plot/SFig5.permutations_HSnCFW.pdf"

## Get data
fullHS = read.csv(fullfileHS, sep ="\t", header = T)

obj = load(permfileHS)
permHSbi = res$VCs; rm(res)
# keep only pheno of interest
fullHS = fullHS[fullHS$taskID %in% permHSbi$taskID[1], ]


## CFW
fullCFW = read.csv(fullfileCFW, sep ="\t", header = T)

obj = load(permfileCFW)
permCFWbi = unique(res$VCs); rm(res)

fullCFW = fullCFW[fullCFW$taskID %in% permCFWbi$taskID[1], ]

## plot 
dict = c("prop_Ad1" = "DGE1", 
         "prop_As1" = "IGE1", 
         "prop_Ad2" = "DGE2", 
         "prop_As2" = "IGE2",
         "corr_Ad1s1" = "cor(DGE1,IGE1)",
         "corr_Ad1d2" = "cor(DGE1,DGE2)",
         "corr_Ad1s2" = "cor(DGE1,IGE2)",
         "corr_Ad2s1" = "cor(DGE2,IGE1)",
         "corr_Ad2s2" = "cor(DGE2,IGE2)",
         "corr_As1s2" = "cor(IGE1,IGE2)"
)

coolors = c("real data"  = adjustcolor("#0038a7", alpha.f = 0.8), #"DGE+IGE+CE"
            "permutations" = "#ededeb")
pch = c("real data"  = 15, #"DGE+IGE+CE"
        "permutations" = 22)

my_boxplot = function(full, perm, vc, pos.lg = "topright", main = "", rotatex= F, line = NULL, leg_cex = 0.8, gap_in = 0.1, ...){
  # --- rotated labels: compute the bottom margin (in inches) ---
  lab = if (is.null(names(vc))) vc else names(vc)
  label_off = 0.1
  lab_h_in = 0.5 # gap between axis and labels
  par(las = 1)
  if(rotatex){
    lab_h_in  = max(strwidth(lab, units = "inches")) * sin(pi / 4) +
      strheight("A", units = "inches") * cos(pi / 4)
    bottom_in = label_off + lab_h_in + 0.5            # + room for the sub text
    
    op = par(las = 1, mai = c(bottom_in, par("mai")[2:4]))
    on.exit(par(op), add = TRUE)                      # restore par on exit
  }

  # 1. Reserve a right margin just wide enough for the legend (in inches)
  leg_w_in = max(strwidth(names(coolors), units = "inches", cex = leg_cex)) +
    par("cin")[1] * leg_cex * 4
  op = par(mai = c(par("mai")[1:3], leg_w_in + gap_in + 0.1))
  on.exit(par(op), add = TRUE)
  
  ylimi = range(full[,vc], perm[,vc])
  boxplot(perm[,vc], ylim = ylimi, xaxt = "n", col = coolors["permutations"], 
          main = main, ylab = "estimated value")

  if(rotatex){
    # ticks only, then rotated labels
    axis(1, at = seq_along(vc), labels = FALSE)
    usr   = par("usr")
    y_lab = grconvertY(grconvertY(usr[3], "user", "inches") - label_off, "inches", "user")
    text(x = seq_along(vc), y = y_lab, labels = lab, srt = 45, adj = c(1, 1), xpd = NA)
  } else{
    axis(1, at = seq_along(vc), labels = lab)
  }
  # sub text below the rotated labels
  in_per_line = par("mai")[1] / par("mar")[1]
  mtext(paste0(nrow(perm), " permutations"),
        side = 1, line = (label_off + lab_h_in) / in_per_line + 1)

  points(x = seq_along(vc), full[,vc], col = coolors["real data"], pch = pch["real data"], cex = 1.5)
  #points(x = seq_along(vc), full[,vc], col = coolors["DGE+IGE+CE"], pch = pch["DGE+IGE+CE"], cex = 1.5)
  
  if(!is.null(line)){
    abline(h = line, lty = 2, lwd = 1.5, col = "#109e0a")
  }#else{
  #  abline(h = 0, lty = 2, lwd = 1.5, col = "grey")
  #}
  
  gap_in = 0.1
  op2 = par(xpd = NA)
  on.exit(par(op2), add = TRUE)
  usr   = par("usr")
  x_leg = grconvertX(grconvertX(usr[2], "user", "inches") + gap_in, "inches", "user")
  
  legend(x = x_leg, y = usr[4], xjust = 0, yjust = 1,
         legend = names(coolors),
         pch    = pch[names(coolors)],
         col    = c(coolors[-length(coolors)], "black"),
         pt.bg  = c(NA, coolors[length(coolors)]),
         pt.lwd = c(NA, 1),
         lwd = c(0, 1.2),
         pt.cex = 1.2,
         cex = leg_cex,
         bty    = "n")
}

pdf(outpdf, h = 5, w = 8)
par(mfrow = c(1,2))
vc = c("prop_As1", "prop_Ad2")
names(vc) = dict[vc]
my_boxplot(full = fullHS, perm = permHSbi, main = "HS mice", 
           vc = vc, pos.lg = "topleft", boxwex = 0.5)
vc = c("corr_Ad2s1")
names(vc) = dict[vc]
my_boxplot(full = fullHS, perm = permHSbi, main = "HS mice", 
           vc = vc, pos.lg = "topleft", boxwex = 0.5)

vc = c("prop_As1", "prop_Ad2")
names(vc) = dict[vc]
my_boxplot(full = fullCFW, perm = permCFWbi,
           vc = vc, main = "CFW mice", 
           boxwex = 0.5)

vc = c("corr_Ad2s1")
names(vc) = dict[vc]
my_boxplot(full = fullCFW, perm = permCFWbi, main = "CFW mice", 
           vc = vc, pos.lg = "topleft", boxwex = 0.5)
dev.off()



