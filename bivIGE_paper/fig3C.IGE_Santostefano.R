# Defining new classes according to our definition
new_class = c(                                "aggression" = "social behaviour",      
                                           "approach rate" = "social behaviour",      
                                  "naso-anal contact rate" = "social behaviour",      
                                           "mounting rate" = "social behaviour",      
                                            "rearing rate" = "social behaviour",      
                             "reciprocal latency to fight" = "social behaviour",       
                                         "neck bite marks" = "assumed\nsocial behaviour",      
                                         "tail bite marks" = "assumed\nsocial behaviour",      
                                         "body bite marks" = "assumed\nsocial behaviour",      
                                           "feeding rates" = "non-social\nbehaviour",      
                        "dyadic fights (fighting ability)" = "social behaviour",      
                                     "old cone hoard size" = "non-social\nbehaviour",      
                                     "new cone hoard size" = "non-social\nbehaviour",  
                        "Performing feather pecks 6 weeks" = "social behaviour",          
               "Performing bouts of feather pecks 6 weeks" = "social behaviour",      
                       "Performing feather pecks 38 weeks" = "social behaviour",
              "Performing bouts of feather pecks 38 weeks" = "social behaviour",      
                       "Performing feather pecks 69 weeks" = "social behaviour",
              "Performing bouts of feather pecks 69 weeks" = "social behaviour",      
                                "Running duration - day 1" = "non-social\nbehaviour", 
                                "Running duration - day 2" = "non-social\nbehaviour",  
                                "Running duration - day 3" = "non-social\nbehaviour",
                                "Running duration - day 4" = "non-social\nbehaviour",  
                                "Running duration - day 5" = "non-social\nbehaviour",
                                 "Running duration- day 6" = "non-social\nbehaviour",  
                                   "Running speed - day 1" = "non-social\nbehaviour",
                                   "Running speed - day 2" = "non-social\nbehaviour",  
                                   "Running speed - day 3" = "non-social\nbehaviour",
                                   "Running speed - day 4" = "non-social\nbehaviour",  
                                    "Running speed- day 5" = "non-social\nbehaviour",
                                    "Running speed- day 6" = "non-social\nbehaviour",
                                 "Feather condition score" = "assumed\nsocial behaviour",
                                        "social dominance" = "social behaviour",
                                 "nesting site preference" = "non-social\nbehaviour",  
                            "neck feather condition score" = "assumed\nsocial behaviour",      
                            "back feather condition score" = "assumed\nsocial behaviour",
                            "rump feather condition score" = "assumed\nsocial behaviour",      
                           "belly feather condition score" = "assumed\nsocial behaviour",
                                                 "divorce" = "social behaviour",      
                                              "Aggression" = "social behaviour",      
                                        "Social dominance" = "social behaviour")
coolors = c("assumed\nsocial behaviour" = "grey65", 
            "social behaviour" = "grey65", 
            "non-social\nbehaviour" = "#B894B1", 
            "non-behaviour" = "#70928D") 

# Loading plot data stored with fig3C.get_data_Santostefano.R
git_dir = "~/git/others/meta-analysis_IGEs/"
load("./plot/rev_Santostefano/meta_model_IGE_subset1A_plot.Rdata") # IGE.subset1A_plot
#plot(y = IGE.subset1A_plot$x, x = IGE.subset1A_plot$y, cex = IGE.subset1A_plot$size)

dataset.IGE.subset1A <- read.csv(file.path(git_dir,"data/subsets/dataset_IGE_subset1A.csv"),
                                 header=T)
data.1A = cbind(dataset.IGE.subset1A[,c("Paper_id","Record_id", "Trait_name", "Trait_category", "Social_h2_2")], IGE.subset1A_plot)
#summary(data.1A$y - data.1A$Social_h2_2)

lev = c("non-behaviour", "non-social\nbehaviour","social behaviour", "assumed\nsocial behaviour")
data.1A[,"new_category"] = factor(sapply(data.1A[,"Trait_name"], function(x){ifelse(x %in% names(new_class), new_class[x], "non-behaviour")}), 
                                  levels = lev)
revlev = lev
data.1A[,"tick"] = as.numeric(factor(data.1A$new_category, levels = revlev))

# Changing dimension 
# precision: true min = 11
#            true max = 182
old_min <- min(data.1A$size)
old_max <- max(data.1A$size)
new_min <- 0.6
new_max <- 2

data.1A$size_transf <- (data.1A$size - old_min) / (old_max - old_min) * (new_max - new_min) + new_min

# Setting up colours and labels
at.dict = seq_along(revlev); names(at.dict) = revlev

# Legend prep
# size mapping
lgd.size = c("11" = min(data.1A$size_transf), "182" = max(data.1A$size_transf))

# Define breaks and their labels
#breaks <- c(20, 100, 200)
breaks <- c(50, 100, 150)
names(breaks) = breaks

# Calculate sizes proportional to the range [0.6, 2]
# Linear interpolation between min and max values
min_val <- min(as.numeric(names(lgd.size)))
max_val <- max(as.numeric(names(lgd.size)))
min_size <- min(lgd.size)
max_size <- max(lgd.size)

# Calculate proportional sizes for each break
sizes <- min_size + (breaks - min_val) / (max_val - min_val) * (max_size - min_size)

at.dict = seq_along(lev); names(at.dict) = lev
#coolors = sapply(at.dict, function(x) ifelse(x == 3 | x == 4, "#B894B1", "grey50"))

pdf("./plot/rev_Santostefano/fig3.IGE_beh_nonb.pdf", h=6, w=5)
par(mar=c(6.1,4.6,3.1,6.6))

set.seed(2)
plot(jitter(data.1A$tick, amount = 0.1), 
     data.1A$Social_h2_2,
     xaxt="n",
     ylab = "IGE",
     xlab = "", main = "Santostefano et al. (2025)",
     xlim = c(min(at.dict) - 0.1, max(at.dict)+0.1),
     ylim = c(0, 0.32), # 0.32 as in our data
     col = coolors[as.character(data.1A$new_category)], pch=21,
     bg = adjustcolor(coolors[as.character(data.1A$new_category)], alpha = 0.5),
     cex = data.1A$size_transf, las=1,
     cex.axis = 1.2, cex.lab = 1.4)

axis(1, at = at.dict, labels = F)
#axis(1, at = at.dict, labels = names(at.dict), las =1,  cex.axis = 1.2, line = 1, lwd=0, tick = F)
text(x = at.dict,
     y = par("usr")[3]-abs(par("usr")[3]/1.5),
     adj = 1,
     labels = names(at.dict),
     xpd = T,
     ## Rotate the labels by 35 degrees.
     srt = 40,
     cex = 1.2)

# Draw the legend
legend(par("usr")[2]-par("usr")[2]/20, par("usr")[4], # position
       legend = names(breaks),         # text labels
       pch = 21,                       # filled circles
       pt.cex = sizes,                 # size of points
       cex = 0.8,
       title = "Precision (1/SE)",     # legend title
       title.font = 2,
       title.adj = 3,
       title.cex = 0.9,
       bty = "n",                      # no box around legend
       xpd = T, 
       inset = c(-0.2,0))
dev.off()
