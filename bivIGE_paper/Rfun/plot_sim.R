### FUNCTIONS PART 1 ####
###### General Plot functions ####
# F(X)a: Function to plot CI as arrows, on top of boxplots - takes same df as the one for which done boxplot
# USED #
plot_ci_df <- function(df, mode,col="#E64B35FF", ...){
  means = sapply(df, mean) # means # mean for estimates' distribution for each col in df - i.e. each param
  
  if (mode == "sd"){
    sd = sapply(df, sd)
    ci_man_right = means + 1.96*(sd) 
    ci_man_left = means - 1.96*(sd) 
  }else if (mode == "STE"){
    sd = sapply(df, sd)
    ci_man_right = means + 1.96*(sd/sqrt(nrow(df)) )
    ci_man_left = means - 1.96*(sd/sqrt(nrow(df)) )
  }else{
    stop("invalid mode, chose between 'sd' or 'STE'")
  }
  
  cis <- cbind(ci_man_left, ci_man_right)
  #cis3 <- cbind(ci_man_left, ci_man_right)
  #cis2 <- t(sapply(df, function(x) t.test(x)$conf.int)) # t.test has alt hypo: true mean is not equal to 0
  # cis2 is the same as cis3 - this should be commented and only cis should be kept
  
  #cis # sides of the CI for each col in df, null hypo: true mean equal to 0 
  n <- ncol(df)
  
  ## Only one of the two following cis is going to stay
  for (i in 1:n) {
    # i = col.number;  means[i] = mean for col[i]; # cis[i,] = CI for col[i];  
    # start_pt = (i, means[i]); end_pt1(i, cis[i, 1] and end_pt2(i, cis[i, 2])
    # angle = angle
    # code = type of arrow - 1 is flat
    # lwd = line width
    tryCatch({arrows(i+0.1, means[i], i+0.1, cis[i, ], angle = 90, code = 2, length = 0.03, col = col, lwd = 1.5, ...)}, #lwd = 3)}, #lwd=3 used for plot per wsh munich
             error=function(e){cat("Can't plot arrows - error message:",conditionMessage(e), "\n")}, 
             warning = function(w){cat("Can't plot arrows - warning message:",conditionMessage(w), "\n")})
    points(i+0.1, means[i], pch=16, col=col, ...)
  }
  #for (i in 1:n) {
  #  tryCatch({arrows(i-0.1, means[i], i-0.1, cis2[i, ], angle = 90, code = 3, length = 0.05, col = "darkorange", lwd = 1.5)}, 
  #           error=function(e){cat("Can't plot arrows - error message:",conditionMessage(e), "\n")})
  #}
}

# F(X)b: Function to plot horizontal line at 0
# USED
zeroline <- function(lty.line=1){
  abline(h = 0, col = "darkgrey", lty=lty.line)
}

# General boxfill color
#boxfill = alpha("lightblue", 0.5)
boxfill = adjustcolor("#4DBBD5FF", alpha = 0.25)


###### F(X)1: Comparing EST vs SIM values #####
#        Also used to compare STE vs sd(est)
# ARGUMENTS:
#        'est_df': df with estimates, as: "trait1" "trait2" "params..."
#        'sim_pars': simulated parameters, named vector, with names(sim_pars) corresponding to colnames(est_df)
### example variables
#est_df = uni12
#sim_pars = c(sim_prop, sim_corr_uni)
# USED
vs_sim <- function(est_df, sim_pars){
  # Creating output df - copying input df rownames and colnames
  estVSsim = data.frame(matrix(nrow = nrow(est_df), ncol = ncol(est_df),
                               dimnames = list(rownames(est_df), colnames(est_df)) ))
  #estVSsim[,"trait1"] = est_df[,"trait1"] 
  #estVSsim[,"trait2"] = est_df[,"trait2"] 
  estVSsim[,traits] = est_df[,traits] 
  estVSsim[1:2,]
  
  for (p in colnames(est_df[,-seq_along(traits)]) ){ # col1 and col2 are "trait1" and "trait2"
    sim_val = sim_pars[p]
    print(paste("sim", p, "=", sim_val))
    estVSsim[,p] =  est_df[,p] - rep(sim_val, nrow(estVSsim))
  }
  # estVSsim[1:2,]; est_df[1:2,]; sim_pars
  return(estVSsim)
}


# PLOT FUNCTIONS
###### F(X)2: boxplot 1 est/STE df vs sim/sd(est) value ####
#           NB: other f(x) used: vs_sim(), plot_ci(), zeroline()
#           'par_sim' = named vector with names of param and sim value
# USED # 
library(scales)
bplot_estSim <- function(est1, x, par_sim, ylimi=NULL, lty.line=2, cex.ylab = 1, dots=F,...){
  #par(pch=20, cex.axis= 0.9) # las = 2 to put axes tick labels perpendicular to ax
  est1VSsim <- vs_sim(est1[,c(traits, names(par_sim))], par_sim) # vs_sim(est_df, sim_pars)
  #est1VSsim[1:2,]
  par(mar=c(6.1,5.1,4.1,2.1))
  if(dots == T){
    dotcol = boxfill
    boxcol = alpha("white", 0.5)
    add=T; outl=F
    boxplot(est1VSsim[,names(par_sim)], ylab = "",
            main = paste0(x, ' analysis'), sub = paste0("n pairs ", nrow(est1VSsim)), 
            col= alpha("white", 0), border=alpha("white", 0), ylim = ylimi, 
            xaxt="n", las = 1, cex.sub=0.9, ...)
    
    at.dict = seq_along(names(par_sim)); names(at.dict) = names(par_sim)
    xpt = sapply(1:length(at.dict), function(i){x = est1VSsim[,names(at.dict)[i]]; x = rep(at.dict[i], length(x))})
    colnames(xpt) = names(at.dict); rownames(xpt) = rownames(est1VSsim)
    set.seed(20); points(x = jitter(xpt, factor=0.5),
                         y = as.matrix(est1VSsim[,names(par_sim)]), #jitter(counts_oi[,1], factor=1), 
                         pch=16, cex=1, 
                         col = dotcol)
    
  }else{boxcol=boxfill; add=F; outl=T}
  
  boxplot(est1VSsim[,names(par_sim)], ylab = "",
          main = paste0(x, ' analysis'), sub = paste0("n pairs ", nrow(est1VSsim)), 
          col= boxcol, ylim = ylimi, 
          xaxt="n", las = 1, cex.sub=0.9, 
          add=add, outline = outl,...)
  if(length(grep("STE", names(par_sim))) > 0){ylabi = "standard error - sd(estimates)"
  }else{ ylabi = "estimates - simulated value"}
  
  title(ylab=ylabi, #(mean±sd)", # add again if plotting with CI
        #font.lab=2, 
        cex.lab = cex.ylab, line=3.5)
  axis(1, at = 1:length(par_sim), labels=F)
  text(x = 1:length(par_sim),
       y = par("usr")[3]-abs(par("usr")[3]/15),
       adj = 1,
       #y = -0.3,
       labels = names(par_sim),
       xpd = T,
       ## Rotate the labels by 35 degrees.
       srt = 35,
       cex = 1)
  axis(3, at = 1:length(par_sim), labels=round(par_sim, 2), line=-1, tick = F, cex.axis= 0.8)
  
  zeroline(lty.line = lty.line)
  # plotting CI with mean - using arrow
  #plot_ci_df(est1VSsim[,names(par_sim)], mode = "sd") # df = est1VSsim[,names(par_sim)]
  #return(p) 
}
### example variables
#est1 = uni12
#par_sim = c(sim_prop, sim_corr_uni)
#x = "uni"
### 
