# fx1. select columns 
select_col = function(mx, colname){
  cat("selecting all columns with '", colname,"'\n")
  mx = mx[, grep(colname, colnames(mx), value=T)]
  return(mx)
}