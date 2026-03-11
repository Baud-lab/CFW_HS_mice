pval_wtest <- function(est, se, constr=0, adj=NULL){
  w = (abs(est) - constr) / se
  pval=pchisq(w^2, df = 1, lower.tail = F)
  # using z statistic, results are the same
  #z <- est / se
  #pval <- 2 * pnorm(-abs(z))
  if(!is.null(adj)){
    pval=p.adjust(pval, method = adj)
  }
  return(pval)
}
