## functions needed for my heatmap #####
suppressMessages(library("corrplot"))

# functions to plot with size as p values -----------------
#    I hacked the function, I still need corrplot package for some colors
#    look at corrplot_size.R function - there are some info there - look at the "size_vector"
#    changed following from here https://github.com/johannes-titz/corrplot/commit/9362f6a7c2fda794b5ef8895b77f0b2ff979092a
#    taken from here https://stackoverflow.com/questions/60410366/r-using-corrplot-to-visualize-two-variables-e-g-correlation-and-p-value-usi

## fx1. transform p
# transform to use pvalues as size of squares in corr plot
transform_p <- function(x, siglev) {
  # Function to transform p-values so that I can use them for size of squares
  y = sapply(x,  function(x) ifelse(x <=siglev, -log10(siglev), -log10(x)))
  #y = (y-min(y))/(max(y)-min(y)) # normalizing to be between 0 and 1
  y = y / (-log10(siglev) ) # normalizing to be between 0 and 1
  #y = log10(x)
  #y <- 0.91 - (0.82) * (1 - exp(-3.82 * x))
  y 
}

## fx2. heatmap with size as p and legend
corrplot_legend= function(mat.toplot, p_mat.toplot, plot.title="", nlg=4, siglev = 0.05, pch.col = "black", pch.cex = 1.2, tl.cex=1, cl.cex=1){ #, maxp = 0.5){
  #arg=(...)
  #p_trans = transform_p(as.numeric(t(p_nom[rowkeep, colkeep])))
  p_trans = transform_p(as.numeric(p_mat.toplot), siglev = max(siglev))
  #p_trans[p_mat.toplot > maxp] = 0
  #as.numeric(p_mat.toplot)
  #p_trans[is.na(p_trans)] = 0
  #warning("transformed_pval resulting as NA are put as 0")
  
  plot(as.numeric(p_mat.toplot), p_trans,
       xlab="pmat", ylab="transformed p", main="square size vs p value")
  
  #plot.new()
  layout(matrix(c(1,0,1,2), 2, 2, byrow = TRUE), 
         widths = c(5,1),
         heights = c(1,1))
  #layout.show(n = 2)
  mar = c(1,1,2,1)
  corrplot_size(mat.toplot, p.mat=p_mat.toplot, size_vector = na.omit(p_trans),
           main=plot.title, # title
           method="square", # method
           insig = 'label_sig', sig.level = siglev, # significance
           pch.cex = pch.cex, pch="*", pch.col = pch.col,  # pch of significance
           col = COL2('RdBu', 200)[200:1], cl.pos = 'b', cl.ratio=0.1, cl.align.text="c", cl.cex=cl.cex, # legend
           tl.col = "grey20", tl.cex = tl.cex, # labels
           na.label = 'square', na.label.col = 'grey80', # handling NAs
           mar=mar)#, tl.offset = 5)
  #m1=siglev[2] #which(p_mat.toplot== max(p_mat.toplot[which(p_mat.toplot <= siglev[2])]) )[1]
  #m2=which(p_mat.toplot== max(p_mat.toplot[which(p_mat.toplot <= 0.2)]) )[1]
  m2=which(p_mat.toplot== max(p_mat.toplot[which(p_mat.toplot <= 0.3)]) )[1]
  m3=which(p_mat.toplot== max(p_mat.toplot[which(p_mat.toplot <= 0.5)]) )[1]
  m4=which(p_mat.toplot== max(p_mat.toplot[which(p_mat.toplot <= 0.9)]) )[1]
  #val1 = as.numeric(p_mat.toplot)
  leg.topl = matrix(c(1, p_trans[c(m2, m3,m4)], c(max(siglev), as.numeric(p_mat.toplot)[c(m2,m3,m4)])), ncol=2)
  rownames(leg.topl) = c(paste0("<=",max(siglev)), round(leg.topl[2:nrow(leg.topl),2], 2))
  print(leg.topl)
  #leg.topl = rbind(leg.topl, matrix(rep(0, nrow(mat.toplot)-nrow(leg.topl)), ncol=1))
  vmar = h*2 #nrow(mat.toplot) #nrow(mat.toplot)/h
  corrplot_size(leg.topl[1:nlg,1,drop=F], 
           main=paste0("p(",adj_toplot,") size"), # title
           method="square", # method
           #insig = 'label_sig', sig.level = siglev, # significance
           #pch.cex = 1.2, pch="*", pch.col = 'black',  # pch of significance
           col = "grey40", cl.pos = 'n', cl.ratio=0, cl.align.text="c", cl.cex=0, # legend
           tl.col = "grey20", tl.cex = 1, tl.srt = 0, tl.offset = 1,tl.pos = "l", # labels
           mar=c(vmar,1,1,0))
}


## fx3. heatmap without p size
corrplot_unsized = function(mat.toplot, plot.title="", methad = "circle", pch.col = "black", pch.cex = 1.2, tl.cex=1, cl.cex=1){ #, maxp = 0.5){
  #arg=(...)
  mar = c(1,1,2,1)
  #sizes = mat.toplot
  #sizes[!is.na(sizes)] = 0.9
  corrplot_size(mat.toplot, #size_vector = na.omit(sizes),
                main=plot.title, # title
                method=methad, # method
                col = COL2('RdBu', 200)[200:1], cl.pos = 'b', cl.ratio=0.1, cl.align.text="c", cl.cex=1, # legend
                tl.col = "grey20", tl.cex = tl.cex, # labels
                na.label = 'square', na.label.col = 'grey80', # handling NAs
                mar=mar) #, tl.offset = 5)
}
