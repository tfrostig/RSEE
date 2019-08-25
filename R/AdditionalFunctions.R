#--------------------------------------
# ADDITIONAL FUNCTIONS FOR SMOOTHING AND SEGMENTATION

# finds homecage coordinates
# returns a vector of length two - c(x.home,y.home)
find.home <- function(x,y,n=5000) {
  s <- sample(x[!is.na(x)],n)
  t <- table(s)
  x.home <- as.numeric(names(t[which.max(t)]))
  y.home <- y[min(which(x==x.home))]
  return(c(x.home,y.home))
}

# fix segments of NAs to be on a straight line between the previous and next point
fix.na <- function(A,home=NULL) {
  x <- A$x
  y <- A$y
  ind.na <- which(is.na(x))
  if (length(ind.na)==0) return(A)
  a.start <- ind.na[which(diff(c(-1,ind.na))>1)]
  a.stop <- ind.na[which(diff(c(ind.na,length(x)+2))>1)]
  n.start <- length(a.start)
  for (i in 1:n.start) {
    ind <- a.start[i]:a.stop[i]
    if (i==1) {
      x[ind] <- x[a.stop[i]+1]
      y[ind] <- y[a.stop[i]+1]
    } else if (i==n.start) {
      x[ind] <- x[a.start[i]-1]
      y[ind] <- y[a.start[i]-1]
    } else {
      if (is.null(home)) {
        x[ind] <- straight.line(x[(a.start[i]-1):(a.stop[i]+1)])[2:(length(ind)+1)]
        y[ind] <- straight.line(y[(a.start[i]-1):(a.stop[i]+1)])[2:(length(ind)+1)]
      } else {
        if (y[a.start[i]-1]==home[2]) { # NAs follow home
          x[ind] <- home[1] ; y[ind] <- home[2]
        } else {
          if (y[a.stop[i]+1]==home[2]) {
            x[ind] <- straight.line(x[c(a.start[i]-1,ind,1)])[2:(length(ind)+1)]
            y[ind] <- straight.line(y[c(a.start[i]-1,ind,1)])[2:(length(ind)+1)]
          } else {
            x[ind] <- straight.line(x[(a.start[i]-1):(a.stop[i]+1)])[2:(length(ind)+1)]
            y[ind] <- straight.line(y[(a.start[i]-1):(a.stop[i]+1)])[2:(length(ind)+1)]
          }
        }
      }
    }

  }
  A$x <- x
  A$y <- y
  return(A)
}

# force all points between the first and last point to a straight line between them
straight.line <- function(x) {
  n <- length(x)-2
  if (n<=0) return (x)
  return (c(x[1],((n:1)*x[1]+(1:n)*x[n+2])/(n+1),x[n+2]))
}
