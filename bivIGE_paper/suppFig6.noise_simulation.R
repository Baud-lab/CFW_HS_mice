### plot simulations with increasing DEE2
sets = c("0.1" = "./data/CFW/simulations/bivariate/set0_DG2_IG1_0.1_s30_estNste.Rdata", 
         "0.44" = "./data/CFW/simulations/bivariate/set0_DEE2_0.44_s20_estNste.Rdata", 
         "2.2" = "./data/CFW/simulations/bivariate/set0_DEE2_2.2_s20_estNste.Rdata"
  
)
simfiles = c("0.1" = "./data/CFW/simulations/bivariate/params_bi_V0.1_S30.txt",
             "0.44" = "./data/CFW/simulations/bivariate/params_bi_V0.44_S20.txt",
             "2.2" = "./data/CFW/simulations/bivariate/params_bi_V2.2_S20.txt"
)

outpdf = "./plot/CFW/SFig6.simulations_noise.pdf"


all_sets = vector("list", length = 3); names(all_sets) = names(sets)
all_simvc = vector("list", length = 3); names(all_simvc) = names(simfiles)
for(n in names(sets)){
  obj = load(sets[n])
  if(n == "0.1"){
    sim = VCs; rm(VCs)
  }else{
    sim = res$VCs; rm(res)
  }
  sim_params = read.csv(simfiles[n], sep = "\t", header = T)
  prop_names = c("var_Ad1","var_As1", "var_Ad2", "var_As2"); corr_names = grep("corr_A", rownames(sim_params), value = T) #c("corr_Ad1d2.rho_Ad1d2")
  vc_sim = c(sim_params[prop_names,"prop_params"], sim_params[corr_names, "set_params"])
  names(vc_sim) = c(gsub("var_", "prop_", prop_names), unlist(sapply(strsplit(corr_names, "\\."), "[[", 1)) )
  all_sets[[n]] = sim
  all_simvc[[n]] = vc_sim
  newn = round(sim_params["var_Ed2","prop_params"], 2)
  names(all_sets)[which(names(all_sets) == n)] = newn
  names(all_simvc)[which(names(all_simvc) == n)] = newn
  rm(sim_params, prop_names, corr_names, vc_sim, sim)
}


# Setting up to plot
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

## Good plotting 
multi_box = function(allest, vc, vc_sim,
                     coolors  = adjustcolor(c("#0038a7", "#724e94", "#d6006f", "#ededeb"), alpha.f = 0.8),
                     sim_cool = "#109e0a", lwd = 2, lty = 2, seglen=0.8,
                     leg_cex = 1, gap_in = 0.1, rotatex = F) {

  # --- rotated labels: compute the bottom margin (in inches) ---
  lab = if (is.null(names(vc))) vc else names(vc)
  label_off = 0.1                                   # gap between axis and labels
  lab_h_in  = max(strwidth(lab, units = "inches")) * sin(pi / 4) +
    strheight("A", units = "inches") * cos(pi / 4)
  bottom_in = label_off + lab_h_in + 0.5            # + room for the sub text
  
  op = par(las = 1, mai = c(bottom_in, par("mai")[2:4]))
  on.exit(par(op), add = TRUE)                      # restore par on exit
  
  # ---- datasets ----
  if (is.null(names(allest))) names(allest) = paste0("data", seq_along(allest))
  cols  = vc
  n_df  = length(allest)
  n_col = length(cols)
  
  # ---- sim values: one named vector, or a list with one per dataset ----
  if (!is.list(vc_sim)) vc_sim = list(vc_sim)
  n_sim = length(vc_sim)
  if (!n_sim %in% c(1, n_df))
    stop("vc_sim must have length 1 or one element per dataset (", n_df, ")")
  sim_vals = lapply(vc_sim, function(v) unname(v[cols]))
  
  # ---- colors: take the first n, or interpolate if there are too few ----
  adapt_cols = function(cl, n) {
    if (length(cl) >= n) cl[seq_len(n)] else colorRampPalette(cl, alpha = TRUE)(n)
  }
  coolors = adapt_cols(coolors, n_df)
  if (!length(sim_cool) %in% c(1, n_sim))
    stop("sim_cool must have length 1 or one color per vc_sim (", n_sim, ")")
  sim_cool = rep_len(sim_cool, n_sim)
  
  # ---- data and positions ----
  # flat list ordered: col1-df1, col1-df2, ..., col2-df1, ...
  data_list = unlist(lapply(cols, function(col) lapply(allest, "[[", col)),
                     recursive = FALSE)
  x_c   = (seq_len(n_col) - 1) * (n_df + 1)      # offset of each column group
  pos   = unlist(lapply(x_c, function(x) x + seq_len(n_df)))
  fills = rep(coolors, times = n_col)
  
  # ---- legend entries ----
  # a single "sim value" entry if all sim lines share a color, otherwise one per dataset
  if (length(unique(sim_cool)) == 1) {
    sim_leg_lab = "simulated value"
    sim_leg_col = sim_cool[1]
  } else {
    sim_leg_lab = paste(names(allest), "sim")
    sim_leg_col = sim_cool
  }
  n_leg_sim  = length(sim_leg_lab)
  leg_labels = c(names(allest), sim_leg_lab)
  
  # 1. Reserve a right margin just wide enough for the legend (in inches)
  leg_w_in = max(strwidth(leg_labels, units = "inches", cex = leg_cex)) +
    par("cin")[1] * leg_cex * 4
  op = par(mai = c(par("mai")[1:3], leg_w_in + gap_in + 0.1))
  on.exit(par(op), add = TRUE)
  
  # 2. Boxplot (ylim includes the sim values so the lines are never clipped)
  boxplot(data_list, at = pos, col = fills, xaxt = "n", xlab = "", ylab = "estimated value",
          ylim = range(c(unlist(data_list), unlist(sim_vals)), na.rm = TRUE))
          #sub  = paste0("n phenotypes: ",
          #              paste0(names(allest), ":", sapply(allest, nrow), collapse = ", ")))
  
  # label each group in the middle
  centers = x_c + (n_df + 1) / 2
  if(rotatex){
    # ticks only, then rotated labels
    axis(1, at = centers, labels = FALSE)
    usr   = par("usr")
    y_lab = grconvertY(grconvertY(usr[3], "user", "inches") - label_off, "inches", "user")
    text(x = centers, y = y_lab, labels = lab, srt = 45, adj = c(1, 1), xpd = NA)
    
  } else{
    axis(1, at = centers, labels = if (is.null(names(cols))) cols else names(cols))
    #axis(1, at = seq_along(vc), labels = lab)
  }
  # sub text below the rotated labels
  in_per_line = par("mai")[1] / par("mar")[1]
  mtext(paste0("n pheno pairs: ", paste0(names(allest), ":", sapply(allest, nrow), collapse = ", ")),
        side = 1, line = (label_off + lab_h_in) / in_per_line + 1, cex = 0.8)
  
  
  # 3. Simulated values
  if (n_sim == 1) {
    # one line across all boxes of each column
    segments(x0 = x_c + 0.5, x1 = x_c + n_df + 0.5,
             y0 = sim_vals[[1]], y1 = sim_vals[[1]],
             lty = lty, lwd = lwd, col = sim_cool)
  } else {
    # one short line per box, matching its dataset
    for (i in seq_len(n_sim)) {
      segments(x0 = x_c + i - 0.4, x1 = x_c + i + 0.4,
               y0 = sim_vals[[i]], y1 = sim_vals[[i]],
               lty = lty, lwd = lwd, col = sim_cool[i])
    }
  }
  
  # 4. Legend outside, aligned with the top-right corner of the plot region
  op2 = par(xpd = NA)
  on.exit(par(op2), add = TRUE)
  usr   = par("usr")
  x_leg = grconvertX(grconvertX(usr[2], "user", "inches") + gap_in, "inches", "user")
  
  legend(x = x_leg, y = usr[4], xjust = 0, yjust = 1,
         legend = leg_labels,
         pch    = c(rep(22, n_df), rep(NA, n_leg_sim)),
         pt.bg  = c(coolors, rep(NA, n_leg_sim)),          # fill of the squares
         col    = c(rep("black", n_df), sim_leg_col),      # square border / line color
         lty    = c(rep(NA, n_df), rep(lty, n_leg_sim)),
         lwd    = c(rep(NA, n_df), rep(lwd, n_leg_sim)),
         seg.len = c(rep(1, n_df), rep(seglen, n_leg_sim)),
         bty    = "n", pt.cex = 2, cex = leg_cex
         )
  
  invisible(NULL)
}

names(all_sets) = paste0("DEE2=",names(all_sets))
names(all_simvc) = paste0("DEE2=",names(all_simvc))

pdf(outpdf, h = 6, w = 7)
vcs = grep("prop",names(all_simvc[[3]]), value = T)
names(vcs) = dict[vcs]
multi_box(allest = all_sets,
          vc_sim = all_simvc,
          vc = vcs, 
          coolors = adjustcolor(c("#ffe500", "#fd8c00", "#fe0000"), alpha.f = 0.5), 
          sim_cool = "#004dff", lwd = 3.5, lty = 5, seglen = 0.8
)

vcs = grep("corr",names(all_simvc[[3]]), value = T)
names(vcs) = dict[vcs]
vcs = na.omit(vcs[match(names(dict), vcs)])
multi_box(allest = all_sets,
          vc_sim = all_simvc,
          vc = vcs, 
          rotatex = T,
          coolors = adjustcolor(c("#ffe500", "#fd8c00", "#fe0000"), alpha.f = 0.5), 
          sim_cool = "#004dff", lwd = 3.5, lty = 5, seglen = 0.8
)

vcname = grep("prop_A", colnames(all_sets[[1]]), value = T)
vcs = gsub("prop_A","STE_A",vcname)
names(vcs) = paste0("se(", dict[vcname], ")")
vcs = na.omit(vcs[match(names(dict), vcname)])
multi_box(allest = all_sets,
          vc_sim = all_simvc,
          vc = vcs, 
          rotatex = T,
          coolors = adjustcolor(c("#ffe500", "#fd8c00", "#fe0000"), alpha.f = 0.5), 
          sim_cool = "#004dff", lwd = 3.5, lty = 5, seglen = 0.8)


vcname = grep("corr_A", colnames(all_sets[[1]]), value = T)
vcs = gsub("corr_A","STE_A",vcname)
names(vcs) = paste0("se(", dict[vcname], ")")
vcs = na.omit(vcs[match(names(dict), vcname)])
multi_box(allest = all_sets,
          vc_sim = all_simvc,
          vc = vcs, 
          rotatex = T,
          coolors = adjustcolor(c("#ffe500", "#fd8c00", "#fe0000"), alpha.f = 0.5), 
          sim_cool = "#004dff", lwd = 3.5, lty = 5, seglen = 0.8)
dev.off()





### quick plotting to see per set 
par(mfrow = c(1,3))
sapply(all_simvc, names)

vcs = grep("prop",names(all_simvc[[1]]), value = T)

ylimi = range(all_sets[[1]][, vcs], all_simvc[[1]][vcs], 
              all_sets[[2]][, vcs], all_simvc[[2]][vcs], 
              all_sets[[3]][, vcs], all_simvc[[3]][vcs])
boxplot(all_sets[[1]][, vcs], ylim = ylimi, ylab = "proportion of phenotypic variance explained", 
        sub = paste0("n simulations = ", nrow(all_sets[[1]])), main = paste0("DEE2 = ", names(all_sets)[1]))
points(x = seq_along(vcs), y = all_simvc[[1]][vcs], col = "red", pch = 16, cex = 2)

vcs = grep("prop",names(all_simvc[[2]]), value = T)
boxplot(all_sets[[2]][, vcs], ylim = ylimi, 
        sub = paste0("n simulations = ", nrow(all_sets[[2]])), main = paste0("DEE2 = ", names(all_sets)[2]))
points(x = seq_along(vcs), y = all_simvc[[2]][vcs], col = "red", pch = 16, cex = 2)

vcs = grep("prop",names(all_simvc[[3]]), value = T)
boxplot(all_sets[[3]][, vcs], ylim = ylimi, 
        sub = paste0("n simulations = ", nrow(all_sets[[3]])), main = paste0("DEE2 = ", names(all_sets)[3]))
points(x = seq_along(vcs), y = all_simvc[[3]][vcs], col = "red", pch = 16, cex = 2)


##
vcs = grep("corr",names(all_simvc[[1]]), value = T)
ylimi = range(all_sets[[1]][, vcs], all_simvc[[1]][vcs], 
              all_sets[[2]][, vcs], all_simvc[[2]][vcs], 
              all_sets[[3]][, vcs], all_simvc[[3]][vcs])
boxplot(all_sets[[1]][, vcs], ylim = ylimi, 
        sub = paste0("n simulations = ", nrow(all_sets[[1]])), main = paste0("DEE2 = ", names(all_sets)[1]))
points(x = seq_along(vcs), y = all_simvc[[1]][vcs], col = "red", pch = 16, cex = 2)

vcs = grep("corr",names(all_simvc[[2]]), value = T)
boxplot(all_sets[[2]][, vcs], ylim = ylimi, 
        sub = paste0("n simulations = ", nrow(all_sets[[2]])), main = paste0("DEE2 = ", names(all_sets)[2]))
points(x = seq_along(vcs), y = all_simvc[[2]][vcs], col = "red", pch = 16, cex = 2)

vcs = grep("corr",names(all_simvc[[3]]), value = T)
boxplot(all_sets[[3]][, vcs], ylim = ylimi, 
        sub = paste0("n simulations = ", nrow(all_sets[[3]])), main = paste0("DEE2 = ", names(all_sets)[3]))
points(x = seq_along(vcs), y = all_simvc[[3]][vcs], col = "red", pch = 16, cex = 2)
