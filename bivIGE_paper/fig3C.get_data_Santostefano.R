library("orchaRd")
####### DGE-IGE corr
git_dir = "~/git/others/meta-analysis_IGEs/"

# Loading data from git repo - https://github.com/ASanchez-Tojar/meta-analysis_IGEs.git
load(file.path(git_dir,"data/models/meta_model_IGE_subset1A.Rdata")) # meta.model.IGE.subset1A
dataset.IGE.subset1A <- read.csv(file.path(git_dir,"data/subsets/dataset_IGE_subset1A.csv"),
                                 header=T)

## prepare data for IGE plot - lines taken from https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/9d73be945e0a3102aeae909a3e826158db627ac3/006_figures.r ####
meta.model.IGE.subset1A_plot <- orchard_plot(mod_results(meta.model.IGE.subset1A, 
                                                         mod = "1",
                                                         group = "Paper_id", 
                                                         data = dataset.IGE.subset1A), 
                                             xlab = "Effect size (r)", 
                                             trunk.size = 2, 
                                             branch.size = 2,
                                             alpha = 0.3,
                                             transfm = "tanh",
                                             fill = T)+
  scale_fill_manual(values="grey") +
  scale_colour_manual(values="grey")+
  theme(legend.direction="horizontal", legend.title = element_text(size =8),
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_blank()) +
  labs( x=expression(Social~h^2))
#summary(IGE.subset1A_plot$y - dataset.IGE.subset1A$Social_h2_2)

# build plot
pg = ggplot_build(meta.model.IGE.subset1A_plot)
# store data from the plot
IGE.subset1A_plot = pg@data[[1]]
# save data to file
save(IGE.subset1A_plot, file="./plot/rev_Santostefano/meta_model_IGE_subset1A_plot.Rdata")
