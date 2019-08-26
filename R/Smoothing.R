#' RSEE.smooth - performs LOWESS and RRM
#' @param object data frame with x and y coordinates.
#'  name of x coordinates vector should contain "x" or "X"
#'   column name of y coordinates vector should contain "y" or "Y"
#' @param path if object is NULL, loads a CSV file of XY coordinates
#' @param project project name
#' @param individual name of individual animal
#' @param forced set to TRUE if no fixed homebase, and FALSE if fixed homebase (free exploration)
#' @param fps frames per second (of the tracked video)
# LOWESS arguments:
#' @param h.loess half window width of LOWESS
#' @param d.loess degree of polynomial (only d=2 implemented for velocities and accelerations)
#' @param r.loess number of LOWESS iterations (recommended: 2)
# RRM arguments:
#' @param frames.rrm RRM arguments (see NOTE)
#' @param cutoff.rrm RRM arguments (see NOTE)
#' @param h.rrm sequence of half-window widths for repeated RM
#' @note  an arrest is a sequence of at least "frames.rrm" frames in which the maximum distance
#         betwee the two farthest points (after RRM) is no more than "cutoff.rrm" cm

# VALUES:
#' @return  The function returns a list with following objects:
#'  \describe{
#'  \item{smooth.data}{ data frame, containing for each frame: smooth x and y coordinates,
#'    smooth velocities and accelerations in x and y, and entry.}
#'  \item{speed.acc}{data frame of smooth point speeds and accelerations}
#'  \item{entry}{data frame of start and end frames for each entry (only 1 entry if forced = T)}
#'  \item{motion.segments}{data frame of start and end frames for each arrest}
#'  \describe{info - list of general information about the session:
#'  \item{project, individual}{input names}
#'  \item{session.length.minutes}{length of the session in minutes}
#'  \item{home}{fixed homebase coordinates (if forced = F)}
#'  \item{number.of.entries}{more than 1 if forced = F}
#'  }
#'  }
#' @export

RSEE.smooth <- function(object=NULL,path=NULL,project = "",individual = "",
                        forced=T,fps=25,h.loess=10,d.loess=2,r.loess=2,
                        frames.rrm=5,cutoff.rrm=1e-6,h.rrm=c(7,5,3,3)) {
  options(warn=-1)
  library(utils)

  # loading track data and fixing NAs
  if (is.null(object)) {
    if (is.null(path)) { stop("No object or path")
    } else object <- read.csv(path)
  }

  which.x <- union(grep("X",names(object)),grep("x",names(object)))
  which.y <- union(grep("Y",names(object)),grep("y",names(object)))
  if ((length(which.x)==0)||(length(which.y)==0)) stop("No X and Y coordinates")

  cat("RSEE path smoother:\n",sep="")

  cat("    Preparing data and Fixing NAs... ",sep="")


  A <- data.frame(x=as.numeric(as.character(object[,which.x])),
                  y=as.numeric(as.character(object[,which.y])))
  rm(object)

  if (!forced) { home <- find.home(A$x,A$y)
  } else home <- NULL
  A <- fix.na(A,home)
  n <- nrow(A)
  x <- rep(0,n)
  y <- rep(0,n)
  v.x <- rep(0,n)
  v.y <- rep(0,n)
  a.x <- rep(0,n)
  a.y <- rep(0,n)

  # find entries
  if (!forced) {
    ind.entry <- which((A$x!=home[1])|(A$y!=home[2]))
    df <- diff(ind.entry)
    ent.start <- ind.entry[which(c(2,df)>1)]
    ent.stop <- ind.entry[which(c(df,2)>1)]
    rm(df)
  } else {
    ind.entry <- 1:n
    ent.start <- 1
    ent.stop <- n
  }
  n.ent <- length(ent.start)
  which.entry <- rep(0,n)
  for (i in 1:n.ent) which.entry[ent.start[i]:ent.stop[i]] <- i

  cat("Done.\n    Running LOESS... ",sep="")
  if (n.ent>1) pb <- txtProgressBar(min = 1, max = n.ent, initial = 1, char = "=",style=3)
  # Performing LOESS separately for each entry
  for (j in 1:length(ent.start)) {
    ind <- ent.start[j]:ent.stop[j]
    temp.x <- LOWESS(A$x[ind],fps,h.loess,d.loess,r.loess)
    x[ind] <- temp.x$x
    v.x[ind] <- temp.x$v
    a.x[ind] <- temp.x$a
    temp.y <- LOWESS(A$y[ind],fps,h.loess,d.loess,r.loess)
    y[ind] <- temp.y$x
    v.y[ind] <- temp.y$v
    a.y[ind] <- temp.y$a
    rm(temp.x,temp.y,ind)
    if (n.ent>1) setTxtProgressBar(pb, j)
  }

  # Performing RRM separately for each entry
  # arrests segments are added
  cat("Done.\n    Running RRM, finding arrests and motion segments... ")
  arr.start <- numeric() # start frames of arrests
  arr.stop <- numeric() # end frames of arrests
  arr.entry <- numeric() # in which entry the arrests was performed

  if (n.ent>1) pb <- txtProgressBar(min = 1, max = n.ent, initial = 1, char = "=",style=3)
  for (j in 1:length(ent.start)) {
    ind <- ent.start[j]:ent.stop[j]
    rrm <- RRM(A$x[ind],A$y[ind],l=frames.rrm,eps=cutoff.rrm,h.seq=h.rrm)
    arr.start <- c(arr.start,ind[rrm$start])
    arr.stop <- c(arr.stop,ind[rrm$stop])
    arr.entry <- c(arr.entry,rep(j,length(rrm$start)))
    rm(rrm,ind)
    if (n.ent>1) setTxtProgressBar(pb, j)
  }

  # set velocity and acceleration to zero in arrests
  # forces X and Y coordinates on a straight line

  for (j in 1:length(arr.start)) {
    ind <- arr.start[j]:arr.stop[j]
    v.x[ind] <- 0
    a.x[ind] <- 0
    x[ind] <- straight.line(x[ind])
    v.y[ind] <- 0
    a.y[ind] <- 0
    y[ind] <- straight.line(y[ind])
  }

  v <- sqrt(v.x^2+v.y^2)
  # find motion segments

  motion.start <- numeric()  # start frames of motion segments
  motion.stop <- numeric() # end frames of motion segments
  motion.entry <- numeric() # which entry the motion segment
  max.v <- numeric()

  n.arr <- length(arr.entry)
  count <- 0
  flag <- 1
  for (j in 1:length(ent.start)) {
    ind <- ent.start[j]:ent.stop[j]
    count <- count + flag
    flag <- 1
    w.m <- c(count:n.arr)[which(arr.entry[count:n.arr]==j)]
    if ((length(w.m)>0)&&(count<=n.arr)) {
      w <- max(w.m)
      if (arr.start[count]>ent.start[j]) {
        motion.start <- c(motion.start,ent.start[j])
        motion.stop <- c(motion.stop,arr.start[count]-1)
        motion.entry <- c(motion.entry,j)
        max.v <- c(max.v,max(v[ent.start[j]:(arr.start[count]-1)]))
      }
      if (w>count) {
        for (k in count:(w-1)) {
          motion.start <- c(motion.start,arr.stop[k]+1)
          motion.stop <- c(motion.stop,arr.start[k+1]-1)
          motion.entry <- c(motion.entry,j)
          max.v <- c(max.v,max(v[(arr.stop[k]+1):(arr.start[k+1]-1)]))
        }
      }
      if (arr.stop[w]<ent.stop[j]) {
        motion.start <- c(motion.start,arr.stop[w]+1)
        motion.stop <- c(motion.stop,ent.stop[j])
        motion.entry <- c(motion.entry,j)
        max.v <- c(max.v,max(v[(arr.stop[w]+1):ent.stop[j]]))
      }
      count <- w
    } else {
      flag <- 0
      motion.start <- c(motion.start,ent.start[j])
      motion.stop <- c(motion.stop,ent.stop[j])
      motion.entry <- c(motion.entry,j)
      max.v <- c(max.v,max(v[ind]))
    }
    rm(ind)
  }


  # Original smoothing file (other are changed later on)
  smooth.data <- data.frame(x=x,v.x=v.x,a.x=a.x,
                            y=y,v.y=v.y,a.y=a.y,entry = which.entry)

  spd.acc <- data.frame(v = v, a = sqrt(a.x^2+a.y^2))


  # entry data
  entry.dat <- data.frame(start=ent.start,end=ent.stop,length=ent.stop-ent.start+1)

  # data frame of arrests data
  arr.data <- data.frame(start=arr.start,end=arr.stop,length=arr.stop-arr.start+1,entry=arr.entry)

  # data frame of motions segements data
  motion.dat <- data.frame(start=motion.start,end=motion.stop,
                           length=motion.stop-motion.start+1,entry=motion.entry,max.v=max.v)

  # info
  info <- list(project = project, individual = individual,
               length.of.session.minutes = n/(fps*60),home = home,number.of.entries = n.ent)

  rm(A,dat,x,y,v.x,v.y,a.x,a.y,v,
     ent.start,ent.stop,ind.entry,
     arr.stop,arr.start,arr.entry,
     motion.start,motion.stop,motion.entry,max.v,
     n.arr,count,home)
  gc()

  cat("Done.\nComplete.",sep="")
  return(list(smooth.data = smooth.data,
              speed.acc = spd.acc,
              entry = entry.dat,
              arrests = arr.data,
              motion.segments  = motion.dat,
              info = info))
}


#' Locally Weighted Scatterplot Smoothing
#' @description Applies LOWESS on select vector as a function of time
#' @param - fps, frame per second
#' @param h - Half window size, the window size is 2 * h  + 1
#' @param d - Degree of polynomial
#' @param r - Number of repeated iterations
#' @param smooth.para - Smoothing parameter, defines
#'
#' @return loc - the smoothed location
#' @return v - velocity obtained by deriving the polynomial at the specific point
#' @return a - accelration, obtained by deriving the polynomial at the specific point twice


LOWESS <- function(x, fps = 25, h = 10, d = 2,  r = 2, smooth.para = 6) {
  n        <- length(x)
  if (h > n) {
    h <- n
    print('Window size cannot be smaller than number of observation, coercing h = n')
  }
  time.vec  <- seq(-h , h , 1) / fps ### For weighting and derivatives
  ### Creating weight
  w  <- (1 - abs( seq(-h, h, 1) / h) ** 3) ** 3 ### Tricube weight (distance on X axis)
  w  <- t(replicate(n, w))  ### Creating weight vector for each observation

  local.test.mat <- NULL
  for (i in 0:d) {
    local.test.mat <- cbind(local.test.mat, (time.vec) ** i)
  }
  ### Creating local regression X vector
  ### Padding with zeros for start
  for (i in 1:h) {
    w[i, seq(1 , h - i + 1, 1)] <- 0
  }
  ### Padding with zeros for end
  for (i in (n - h + 1):n) {
    w[i, seq(h + (n - i) + 2, 2 * h + 1)] <- 0
  }
  ### Creating residual weighting
  delta <- matrix(1, n, 2 * h + 1) # weights due to distance from the fitted curve
  ### Padding X with zeros
  padx  <- c(rep(0, h), x, rep(0, h))
  #### Iterating
  for (i in 1:r) {
    if (i < r) {
      delta <- iterLowess(weightX = w, weightY = delta, modelmat = local.test.mat, x = padx, h = h, smoothpara = smooth.para)
    }
    if (i == r) { ### Last iteration (return movment properties)
      mov.mat <- lastIterLowess(weightX = w, weightY = delta, modelmat = local.test.mat, x = padx, h = h)
    }
  }
  colnames(mov.mat)[1:(d + 1)] <- c('loc', 'v', 'a', rep(NA, 3 - (d + 1)))[1:(d + 1)]
  if (d > 1) {
    mov.mat[ ,3] <- mov.mat[ ,3] * 2 ## Accelration
  }
  return(mov.mat)
}

#' Repeated running medians
#'
#' @param x - Raw x tracked coordinates
#' @param y - Raw y tracked coordinates
#' @param l - Minimum number of frames to be considered as arrest
#' @param eps - Cutoff value for determining arrests (maximum distance allowed to travel during arrest)
#'
#' @return Table of specfying beginning and ends of arrests
#'
RRM <- function(x,y,l=5,eps=1e-6,h.seq=c(7,5,3,3)) {
  n <- length(x)
  if (n<l) return(data.frame(start=numeric(),stop=numeric()))
  if (max(h.seq)>n) h.seq[which(h.seq>n)] <- n
  x.hat <- x
  y.hat <- y
  for (h in h.seq) {
    x.hat <- runmed(x.hat, 2 * h + 1, endrule = 'constant')
    y.hat <- runmed(y.hat, 2 * h + 1, endrule = 'constant')
  }
  x <- x.hat
  y <- y.hat

  arrests.start <- numeric()
  arrests.stop <- numeric()

  i <- 1
  while (i < n) {
    k <- i+1
    flag <- T
    while (flag) {
      df <- sqrt((x[k]-x[i])^2+(y[k]-y[i])^2)
      if (df>eps) {
        flag <- F
        k <- k-1
      } else {
        if (k == n) {
          flag <- F
        }
        else k <- k+1
      }
    }
    if (k- i + 1 > l) {
      arrests.start <- c(arrests.start,i)
      arrests.stop <- c(arrests.stop,k)
    }
    i <- k + 2
  }

  #df <- rep(0,n)
  #if (l>1) {
  #  for (i in 2:l) {
  #    df <- cbind(df,c(sqrt((x[i:n]-x[1:(n-i+1)])^2+(y[i:n]-y[1:(n-i+1)])^2),rep(0,i-1)))
  #  }
  #}
  return(data.frame('start' = arrests.start,'stop' = arrests.stop)) ## Might cause error, if so, switch to list
}

