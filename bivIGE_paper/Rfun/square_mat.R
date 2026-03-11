## function to get a square matrix from a rectangular one; or from a matrix that has different rows and cols - filling the missing ones with NAs

square.mat = function(mat){
  mat = as.matrix(mat)
  
  # Filling columns
  if(ncol(mat) == 0){
    cat("NB: starting matrix has no columns\n")
  }else if (!any(colnames(mat) %in% rownames(mat))){
    cat("NB: not even one colnames in rownames\n")
  }
  
  for(i in rownames(mat)){
    if(!(i %in% colnames(mat))){
      coli = rep(NA, nrow(mat))
      mat = cbind(mat, coli)
      colnames(mat)[which(colnames(mat) == "coli")] = i
    }
  }
  
  # Filling rows
  if(nrow(mat) == 0){
    cat("NB: starting matrix has no rows\n")
  }else if (!any(rownames(mat) %in% colnames(mat))){
    cat("NB: not even one rownames in colnames\n")
  }
  
  for(i in colnames(mat)){
    if(!(i %in% rownames(mat))){
      rowi = rep(NA, ncol(mat))
      mat = rbind(mat, rowi)
      rownames(mat)[which(rownames(mat) == "rowi")] = i
    }
  }
  
  mat = mat[rownames(mat), rownames(mat)]
  return(mat)
}

# example of mat based on cor_mat
#mat = cor_mat[-sample(1:nrow(cor_mat), 10), -sample(1:nrow(cor_mat), 10)]
#newMat = square.mat(mat)
#rownames(newMat) == colnames(newMat)
