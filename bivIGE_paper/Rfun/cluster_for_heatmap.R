## C. functions for clustering -----------------
cluster.mat = function(mat, distance, method){
  d = dist(mat, method = distance)
  tree = hclust(d, method = method)
  return(tree)
}

plot.tree = function(tree, k, subtitle, colors, cex=0.5){
  # For columns
  #tree = tree_col
  #k = 10
  cl_members <- cutree(tree = tree, k = k)
  pars <- par() # store original pars
  # plot dendrogram
  plot(x = tree, labels =  row.names(tree), cex = cex, xlab="", sub= subtitle)
  #Generate borders around each group
  par(lwd=2, mar=c(0,0,0,0)) # pars to make the line width of rectangles bigger
  rect.hclust(tree = tree, k = k, which = 1:k, border = colors[1:k], cluster = cl_members)
  # reset par
  par(lwd=pars$lwd, mar=pars$mar)
}