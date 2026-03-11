library("here")
sourcefun = "./bivIGE_paper/Rfun" 
source(here(sourcefun, "plot_sim.R"))

# To plot HS mice use these lines 
#opt=list(infile = "~/nf_PRJs/nf-CoreQuantGen/simulations/output/simHSmice/toPlot/toPlot_bi_0.1_s20.Rdata",
#         outdir = "./plot/HSmice/"
#         )

# To plot CFW mice use these lines 
opt=list(infile="~/nf_PRJs/nf-CoreQuantGen/simulations/output/simCFW/toPlot/toPlot_bi_0.1_s30_2406.Rdata",
         outdir= "./plot/CFW/"
)
ylimi = c(-0.35,0.35) # y-limits for estimates plot
se_ylimi = c(-0.035, 0.03) # y-limits for standard errors plot

toplot_file = opt$infile
out_dir = opt$outdir
grip = "A"
eff = "GEN"

#### 2. Loading objects ------
# Loading list of objects saved from 4.simRes_toPlot.R and storing in respective variables
load(toplot_file)

# Loading names 
prop_names = grep(grip, toPlot$prop_names, value=T ) # toPlot$prop_names # NB: they have a different order in uni and bi
corr_names.uni = grep(grip, toPlot$corr_uni_names, value=T ) # old sim
corr_names.bi = grep(grip, toPlot$corr_bi_names, value=T ) # old sim
totv_names = toPlot$totv_names # NB: they have a different order in uni and bi
ste_dict = toPlot$ste_dict
corP_names = toPlot$corP_names

# Now the simulated values - NB: sim_corr specific to bivariate are not present in uni
sim_prop = toPlot$sim_prop
sim_corr = toPlot$sim_corr # These are different in uni and bi! 
sim_totv = toPlot$sim_totv

# Now the dataframes 
stopifnot(all(rownames(toPlot$est12) == rownames(toPlot$est_bi_corrs))); est = cbind(toPlot$est12, toPlot$est_bi_corrs) # old sim
stopifnot(all(rownames(toPlot$STE12) == rownames(toPlot$STE_bi_corrs))); STE = cbind(toPlot$STE12, toPlot$STE_bi_corrs) # old sim
corP = toPlot$corP

smp_sizes = c("sample_size1", "sample_size1_cm", "sample_size2", "sample_size2_cm")
traits = c("trait1", "trait2")


############# PREPARING ESTIMATES ################
#### 1. Selecting only genetic parameters: -----
cat("Params that will be analysed and plotted are: \n")
prop_names; corr_names.uni; corr_names.bi
corr_names = c(corr_names.uni, corr_names.bi)


#### 2. Filtering the datasets and the sim values for parameters want to plot -----
# Datasets:
est = est[,c(traits, smp_sizes, prop_names, corr_names)]
cat("Datasets filtered for params to plot: \n")
cat("'est': \n"); est[1:2,];

# Sim_values:
sim_prop = sim_prop[prop_names]
sim_corr = sim_corr[corr_names] # uni in uni; uni + bi in bi

# Check again everything corresponds 
# and correspond to simulated params -> sim_prop, sim_corr, sim_tot
stopifnot(names(sim_prop) %in% colnames(est))
stopifnot(names(sim_corr) %in% colnames(est))

cat("Sim values filtered for params to plot: \n")
cat("'sim_prop'\n"); sim_prop
cat("'sim_corr'\n"); sim_corr


#### 3. Selecting STE_names for params want to plot -------
#    dict is: names(dict) = STE_names ; dict = par_names
dict_prop = na.omit(ste_dict[match(prop_names, ste_dict)]) #ste_dict[which(ste_dict %in% prop_names)] 
dict_corr = ste_dict[match(corr_names, ste_dict)]

# Filtering STE datasets
STE = STE[, c(traits, smp_sizes, names(dict_prop), names(dict_corr))]; 

# Getting parametric STE as sd(estimated param)
chrcol = which(colnames(est) %in% c(traits, smp_sizes))
empSTE = apply(est[,-chrcol], MARGIN = 2, function(x) sd(x))
names(empSTE) = names( ste_dict[match(names(empSTE), ste_dict)]) #names( ste_dict[ste_dict %in% names(empSTE)] )
empSTE = empSTE[colnames(STE[,-chrcol])]

# Plot estimates
h = 6; w = 7 
out_file=file.path(out_dir, "Sfig1.sim_est.pdf")
cat("Saving plot to ", out_file, "\n")
pdf(out_file, height = h, width = w, bg="white") 
par(pch=16, cex.axis= 1.1)

boxfill= adjustcolor("#4DBBD5FF", alpha.f = 0.4)
bplot_estSim(est1 = est, 
             x = "bi", 
             par_sim = c(sim_prop, sim_corr), 
             ylimi = ylimi, # comment this line if don't want to set y limits
             boxwex = 0.6, cex.ylab=1.3, cex.axis=1.1, cex=0.7, pch=16, dots=T)

nemo=dev.off()

# Plot standard errors
out_file=file.path(out_dir, "Sfig1.sim_se.pdf")
cat("Saving plot to ", out_file, "\n")
pdf(out_file, height = h, width = w)  
par(pch=16, cex.axis= 1.1)

boxfill= alpha("#00A087FF", 0.4)
bplot_estSim(est1 = STE, 
             x = "bi", 
             par_sim = empSTE, 
             boxwex = 0.6, cex.ylab=1.3, cex.axis=1.1, cex=0.7, pch=16,
             ylimi = se_ylimi, # comment this line if don't want to set y limits
             dots=T)  
nemo=dev.off()


# Plot coverage probability
library("rsimsum")
out_file=file.path(out_dir, "Sfig1.sim_covProb.pdf")
pdf(out_file, h = 7, w = 12) 
par(mar = c(4, 8, 3, 2), mfrow=c(1,2))  # Margins: (bottom, left, top, right)
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
    if(length(s) != 0){
      este = merge(est[,p, drop=F], STE[,s, drop=F], by="row.names")
      sme = simsum(este, p, se=s, true=sim_prop[p])
      x = summary(sme)
      rownames(x$summ) = x$summ$stat
      est_ci = x$summ[perf,c("est","lower","upper")]
      return(est_ci)
    }else{
      return(rep(NA,3))
    }
  })
  
  perf_all = as.data.frame(rbind(t(perf_prop),t(perf_corr)))
  perf_all = perf_all[!apply(perf_all, 1, function(x) all(is.na(x))), ]
  
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
       #main = "Bias", font.main = 2)
       main = main, font.main = 2)
  
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

