col.grp <- function(n,b) {
  colorvec <- vector(mode="character", length=length(n))
  plotcols <- rainbow(length(b))
  for (i in 1:length(n)) {
    for (j in 1:length(b)) {
      if ( n[i] == b[j] ) {
        colorvec[i] = plotcols[j]
      }
    }
  }
  c(colorvec)
}
bg.colors <- function(n) {
  colorvec <- vector(mode="character", length=nrow(n))
  r <- rainbow(4)
  for (i in 1:nrow(n)) {
    if (n[i,1] >= 2) {
      if (n[i,2] >= 2) {
        colorvec[i] = "purple"
      }
      else {
        colorvec[i] = "cyan"
        
      }
    }
    if (n[i,1] <= -2) {
      if (n[i,2] >= 2) {
        colorvec[i] = "red"
      }
      else {
        colorvec[i] = "salmon"
      }
    }
    if (n[i,2] >= -2 && n[i,2] <= 2) {
      colorvec[i] = "grey"
    }
    if (n[i,1] >= -2 && n[i,1] <= 2) {
      colorvec[i] = "grey"
    }
}
  c(colorvec)
}

