library("corrplot")
# Plotting only the color legend
# Define the range of values
zlim <- c(-1, 1)

# Create a color palette
color_palette <- colorRampPalette(COL2('RdBu', 200)[200:1])

# Create dummy data for the legend
legend_values <- seq(zlim[1], zlim[2], length.out = 200)

#### Plot HORIZONTAL legend
# aspect ratio: 494, 149
pdf("./plot/fig5.legend_corr_hz.pdf", w = 5, h = 1.5)
par(mar = c(4, 2, 2, 2)) # Adjust margins to make space for the horizontal legend
image(
  x = legend_values, y = 1, z = t(matrix(legend_values, nrow = 1)),
  col = color_palette(200),
  axes = FALSE, xlab = "", ylab = ""
)
box(lwd=2)

## add labels
axis(1, at = seq(-1, 1, by = 0.5), labels = seq(-1, 1, by = 0.5), cex.axis = 2, 
     lwd.ticks = 0, lwd = 0, line=0.5)
## add the ticks
axis(1, at = seq(-1, 1, by = 0.5), labels = rep("", length(seq(-1, 1, by = 0.5))), cex.axis = 2, 
     lwd.ticks = 2, lwd = 2)
dev.off()

## ### Plot VERTICAL legend
## par(mar = c(2, 6, 2, 2)) # Adjust margins to make space for the horizontal legend
## image(
##   y = legend_values, 
##   x = 1, 
##   z = matrix(legend_values, nrow = 1),
##   col = color_palette(200),
##   axes = FALSE, xlab = "", ylab = ""
## )
## box(lwd=2)
## 
## ## add labels
## axis(2, at = seq(-1, 1, by = 0.5), labels = seq(-1, 1, by = 0.5), cex.axis = 2, 
##      lwd.ticks = 0, lwd = 0, line=0.5, las=1)
## ## add the ticks
## axis(2, at = seq(-1, 1, by = 0.5), labels = rep("", length(seq(-1, 1, by = 0.5))), cex.axis = 2, 
##      lwd.ticks = 2, lwd = 2)


# Plotting only square's legend
library("here")
sourcefun = "./bivIGE_paper/Rfun/" # TODO: might have to change this path relative to where code is
source(here(sourcefun, "corrplot_size.R")) # `corrplot_size`: function to get heatmap - based on corrplot and with option "size_vector"

#leg.topl when full DGE-IGE matrix, fdr 0.1
#           [,1]      [,2]
# ≤0.1 1.0000000 0.1000000
# 0.3  0.5254906 0.2982012
# 0.5  0.3018395 0.4990689
# 0.9  0.0459466 0.8996082
leg.topl = matrix(c(1.0000000,0.5254906,0.3018395,0.0459466, 0.1000000,0.2982012,0.4990689,0.8996082), ncol=2)
rownames(leg.topl) = c(" ≤0.1 ", " 0.3 ", " 0.5 ", " 0.9 ")
colnames(leg.topl) = c("", "")
nlg = 4
adj_toplot = "fdr"
print(leg.topl)

pdf("./plot/fig5.legend_FDR_vt.pdf", w=8, h=5)
# Vertical
corrplot_size(leg.topl[1:nlg,1,drop=F], 
              main=paste0("p(",adj_toplot,") size"), cex.main = 3,# title
              method="square", # method
              col = "grey40", cl.pos = 'n', cl.ratio=0, cl.align.text="c", cl.cex=0, # legend
              tl.col = "black", tl.cex = 3, tl.srt = 0, tl.offset = 1, tl.pos = "l", # labels
              mar=c(4,4,4,4))
dev.off()

# Horizontal
#corrplot_size(t(leg.topl)[1,1:nlg,drop=F], 
#              main=paste0("p(",adj_toplot,") size"), cex.main = 3,# title
#              method="square", # method
#              col = "grey40", cl.pos = 'n', cl.ratio=0, cl.align.text="c", cl.cex=0, # legend
#              tl.col = "black", tl.cex = 3, tl.srt = 0, tl.offset = 1, tl.pos = "t", # labels
#              mar=c(4,7,4,7))