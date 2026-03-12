## plot results of simulations for sexvariate analysis
suppressMessages(library("here"))
library("rsimsum") # to plot coverage probability (and bias)

sourcefun = "./bivIGE_paper/Rfun" 
source(here(sourcefun, "plot_sim.R")) 

# To plot CFW results
inrdata = "~/nf_PRJs/nf-CoreQuantGen/simulations/output/simCFW/VD/bivariate/mfphenos/pruned_dosages_DGE_IGE_IEE_cageEffect/pheno_allm_estNste.Rdata"
simP_file = "~/nf_PRJs/nf-CoreQuantGen/simulations/output/simCFW/mockphenos/cage_623/set0/DG2_IG1/0.1/30/params_bi_V0.1_S30.txt"
outfile = "./plot/CFW/SFig3.sim_sexvar_est_se_covProb.pdf"

# To plot HSmice results
inrdata = "~/nf_PRJs/nf-CoreQuantGen/simulations/output/simHSmice/VD/bivariate/mfphenos/Andres_kinship_DGE_IGE_IEE_cageEffect/pheno_allm_estNste.Rdata"
simP_file = "~/nf_PRJs/nf-CoreQuantGen/simulations/output/simHSmice/mockphenos/cage_549/set0/DG2_IG1/0.1/30/params_bi_V0.1_S30.txt"
outfile = "./plot/HSmice/SFig3.sim_sexvar_est_se_covProb.pdf"

ylimi = c(-0.85, 0.59) # this is to have the same axis for the two populations

load(inrdata)
VCs = res$VCs
dim(VCs)

prop_names = grep("A", grep("prop", colnames(VCs), value =T), value=T)
corr_names = grep("A", grep("corr", colnames(VCs), value =T), value=T)

traits = c("trait1", "trait2")
smp_sizes = c("sample_size1", "sample_size1_cm", "sample_size2", "sample_size2_cm")
est = VCs[, c(traits, prop_names, corr_names)]

ste_dict = c(prop_names, corr_names)
names(ste_dict) = c(gsub("prop", "STE", prop_names), gsub("corr", "STE", corr_names) )
STE = VCs[, c(traits, names(ste_dict))]

# Load simulated values
simP = read.delim(simP_file)
sim_prop = simP[gsub("prop", "var", prop_names), "prop_params"]; names(sim_prop) = prop_names
sim_corr = simP[sapply(corr_names, function(x) paste0(x, ".", gsub("corr", "rho", x))), "set_params"]; names(sim_corr) = corr_names

# Plotting boxplot of estimates' differences

pdf(outfile, h=7, w=9)
par(pch=16, cex.axis= 1.1) 
boxfill= alpha("#4DBBD5FF", 0.4)
bplot_estSim(est1 = est, 
             x = paste0("sex-diff (", unique(VCs[,"sample_size1"]), " males; ", unique(VCs[,"sample_size2"]), " females)"), 
             par_sim = c(sim_prop, sim_corr), boxwex = 0.6, cex.ylab=1.3, cex.axis=1.1, cex=0.7, pch=16,
             ylim=ylimi, # comment this line if don't want set y-limits
             dots=T) 

# Plotting boxplot of STEs' differences
chrcol = which(colnames(est) %in% c(traits))
empSTE = apply(est[,-chrcol], MARGIN = 2, function(x) sd(x))
names(empSTE) = names(ste_dict[match(names(empSTE), ste_dict)]) #names( ste_dict[ste_dict %in% names(empSTE)] )
empSTE = empSTE[colnames(STE[,-chrcol])]

boxfill= alpha("#00A087FF", 0.4)
ylimi = c(-0.06, 0.13) # 0.13: removes 3 outliers in HSmice and 1 in CFW 
bplot_estSim(est1 = STE, 
             x = paste0("sex-diff (", unique(VCs[,"sample_size1"]), " males; ", unique(VCs[,"sample_size2"]), " females)"), 
             par_sim = empSTE, 
             boxwex = 0.6, cex.ylab=1.3, cex.axis=1.1, cex=0.7, pch=16,
             dots=T, ylimi = ylimi)  


# bias and coverage prob
par(mar = c(6, 8, 3, 2), mfrow=c(1,2), cex.axis=1)  # Margins: (bottom, left, top, right)

plot_performance = function(perf, main){
  perf_corr <- sapply(names(sim_corr), function(p){
    s=names(ste_dict)[which(ste_dict == p)]
    este = merge(est[,p, drop=F], STE[,s, drop=F], by="row.names")
    sme = simsum(este, p, se=s, true=sim_corr[p])
    x = summary(sme)
    rownames(x$summ) = x$summ$stat
    est_ci = x$summ[perf,c("est","lower","upper")]
    return(est_ci)
  })
  
  perf_prop <- sapply(names(sim_prop), function(p){
    s=names(ste_dict)[which(ste_dict == p)]
    este = merge(est[,p, drop=F], STE[,s, drop=F], by="row.names")
    sme = simsum(este, p, se=s, true=sim_prop[p])
    x = summary(sme)
    rownames(x$summ) = x$summ$stat
    est_ci = x$summ[perf,c("est","lower","upper")]
    return(est_ci)
  })
  
  perf_all = as.data.frame(rbind(t(perf_prop),t(perf_corr)))
  
  variables <- rownames(perf_all)
  estimates <- unlist(perf_all[,"est"])
  lower_ci <- unlist(perf_all[,"lower"])
  upper_ci <- unlist(perf_all[,"upper"])
  if(perf=="bias"){ # to put parentheses at 5% of true value instead of montecarlo SE
    lower_5p = 0 - 0.05 * c(sim_prop, sim_corr)[variables] # or 0.01
    upper_5p = 0 + 0.05 * c(sim_prop, sim_corr)[variables] # or 0.01
    xlimi = range(c(lower_ci, upper_ci, lower_5p, upper_5p))
  }else{
    xlimi = range(c(lower_ci, upper_ci))
  }
  
  # Set up plotting parameters
  y_positions <- length(variables):1  # Reverse order for variables
  
  plot(NULL, xlim = xlimi,
       ylim = c(0.5, length(y_positions) + 0.5),
       xlab = "", ylab = "", yaxt = "n", 
       main = main, font.main = 2, sub=paste0(nrow(est), " sim pheno; ", unique(VCs[,"sample_size1"]), " males; ", unique(VCs[,"sample_size2"]), " females")
  )
  
  # Add y-axis labels
  lobels = c(sapply(names(sim_prop), function(x) paste0(x, " (",round(sim_prop[x], 2), ")")),
             sapply(names(sim_corr), function(x) paste0(x, " (",round(sim_corr[x], 2), ")")))
  axis(2, at = y_positions, labels = lobels[variables], las = 1, cex.axis = 0.9)
  
  # Add vertical reference line at x = 0 - for bias
  if (perf=="bias"){
    abline(v = 0, col = "gray80")
  }else if(perf =="cover"){
    abline(v = 0.95, col = "gray80")
  }
  
  # Add confidence intervals and lollipop elements
  for (i in 1:length(y_positions)) {
    # Horizontal line for confidence interval (thinner, gray line)
    #lines(c(lower_ci[i], upper_ci[i]), c(y_positions[i], y_positions[i]), 
    #      col = "gray50", lwd = 0.8)
    
    # Parentheses for confidence interval
    text(lower_ci[i], y_positions[i], "[", col = "gray50", pos = 2, offset=0)
    text(upper_ci[i], y_positions[i], "]", col = "gray50", pos = 4, offset=0)
    
    # Parentheses for 5% of estimate
    #if(perf=="bias"){
    #  text(lower_5p[i], y_positions[i], "(", col = "darkred", pos = 2, offset=0)
    #  text(upper_5p[i], y_positions[i], ")", col = "darkred", pos = 4, offset=0)
    #  
    #}
    
    # Lollipop stem: line from 0 or 0.95 to the estimate
    if (perf=="bias"){
      lines(c(0, estimates[i]), c(y_positions[i], y_positions[i]), 
            col = "black", lwd = 1.5)  
    }else if(perf =="cover"){
      lines(c(0.95, estimates[i]), c(y_positions[i], y_positions[i]), 
            col = "black", lwd = 1.5)
    }
    
    # Lollipop point
    points(estimates[i], y_positions[i], pch = 16, cex = 1.2)
    
    # Add estimate value as text
    text(estimates[i], y_positions[i] + 0.25, sprintf("%.3f", estimates[i]), 
         cex = 0.8)
  }
}

plot_performance("cover", "coverage probability")
#plot_performance("bias", "bias")
dev.off()

