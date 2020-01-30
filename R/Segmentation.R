#' RSEEsegment - performs segmentation to lingering and progression
#' @param object list of results from RSEEsmooth
#' @param initial.n.gauss number of initial Gaussians to fit (recommended 2)
#' @param which.intersect if zero, returns the intersection between the leftmost and the rightmost Gaussians.
#'   In case of two Gaussians, just the intersection between them)
#'   In case of more than two Gaussians, the user can choose which one (from left to right)
#' @param min.length.motion.seg minimum length (in frames) of motions segments to use for segmentation
#' @param user.permission if TRUE the user gets to verify the segmentation results,
#    and choose number of Gaussians and threshold manuall
#' @return object - updates input list "object" with the following components:
#' \describe{
#' \item{progression.lingering}{an additional data frame including start and end time of each segment,
#'     prog - 1 if progression, 0 if lingering. entry - in which entry was the segment performed}
#' \item{motion.segments}{updated to include classification to lingering and progression (prog = 0 or 1)}
#' \item{info}{updated to include velocity.threshold.cm.s - the velocity threshold between lingering and
#'     progression (cm/s, non transformed)}}
#' @export
RSEEsegment <- function(object,initial.n.gauss=2,which.intersect=0,
                         min.length.motion.seg=2,user.permission=F) {
  options(warn=-1)
  library(utils)
  # calculate gaussian intersection

  motion <- object$motion.segments

  ind.segs <- which(motion$length>=min.length.motion.seg)
  v <- motion$max.v[ind.segs]
  em <- EMgaussian(log(v+2),num.gauss=initial.n.gauss)
  x.intersect <- calcIntersection(log(v+2),obj=em,which.intersect=which.intersect,
                                   save.plot=NULL,x.lab="Maximum segment velocity (log cm/s +2)")

  if (user.permission) {
    n.gauss <- 2
    flag <- T
    while (flag) {
      a <- readline("Is this segmentation ok [Y/N]? ")
      if ((a!="N")&&(a!="n")) {
        flag <- F
      } else {
        cat("Select one of the following:\n   [1] Select intersection manually.",
            "\n   [2] Add Gaussian.\n   [3] Subtract Gaussian.\n   [4] return to original settings.",
            "\n   [other] Keep current settings.",sep="")
        a <- readline(paste("[1] Select intersection manually.", "[2] Add Gaussian.",
                            "  [other] Keep current settings.   "))
        if (a=="1") x.intersect <- locator(1,pch=3,col="dark green",lwd=3,type="p")$x
        if (a=="2") {
          n.gauss <- n.gauss+1
          em <- EMgaussian(log(v+2),num.gauss=n.gauss)
          x.intersect <- calcIntersection(log(v+2),obj=em,which.intersect=which.intersect,
                                           save.plot=NULL,x.lab="Maximum segment velocity (log cm/s +2)")
        }
        if (a=="3") {
          n.gauss <- n.gauss-1
          em <- EMgaussian(log(v+2),num.gauss=n.gauss)
          x.intersect <- calcIntersection(log(v+2),obj=em,which.intersect=which.intersect,
                                           save.plot=NULL,x.lab="Maximum segment velocity (log cm/s +2)")
        }
        if (a=="4") {
          n.gauss <- 2
          em <- EMgaussian(log(v+2),num.gauss=n.gauss)
          x.intersect <- calcIntersection(log(v+2),obj=em,which.intersect=which.intersect,
                                           save.plot=NULL,x.lab="Maximum segment velocity (log cm/s +2)")
        }
      }
    }
  }

  # create segmentation table

  thres <- x.intersect
  thres.non.transformed <- exp(thres)-2
  which.prog <- ind.segs[which(v>thres.non.transformed)]
  motion$prog <- rep(0,nrow(motion))
  motion$prog[which.prog] <- 1
  object$motion.segments <- motion


  prog.ling <- data.frame(start=0,end=0,length=0,prog=0,entry=0)

  entry <- object$entry
  n.ent <- nrow(entry)
  n.prog <- length(which.prog)
  ind <- rep(-1,entry$end[n.ent]) # -1 : no entry. 0: lingering. 1: progression
  ent <- rep(0,entry$end[n.ent])
  for (i in 1:n.ent) {
    ind[entry$start[i]:entry$end[i]] <- 0
    ent[entry$start[i]:entry$end[i]] <- i
  }
  for (i in 1:n.prog) ind[motion$start[which.prog[i]]:motion$end[which.prog[i]]] <- 1

  w.df <- which(c(1,diff(ind))!=0)
  for (i in 1:(length(w.df)-1)) {
    if (ind[w.df[i]]!=(-1)) prog.ling <- rbind(prog.ling,c(w.df[i],w.df[i+1]-1,
                                                           w.df[i+1]-w.df[i],
                                                           ind[w.df[i]],ent[w.df[i]]))
  }
  i <- i+1
  if (ind[w.df[i]]!=(-1)) prog.ling <- rbind(prog.ling,c(w.df[i],entry$end[n.ent],
                                                         entry$end[n.ent]-w.df[i]+1,
                                                         ind[w.df[i]],ent[w.df[i]]))
  prog.ling <- prog.ling[-1,]
  rownames(prog.ling) <- NULL

  object$progression.lingering <- prog.ling

  object$info$velocity.threshold.cm.s <- thres.non.transformed
  cat("Done.")
  return(object)
}





# fits num.gauss gaussians to a vector y
EMgaussian <- function(y,num.gauss=2,mu=NULL,sigma=NULL,p=NULL,max.iter=200,err=0.001) {
  count <- 0
  flag <- T
  n <- length(y)
  if (is.null(p)) p <- matrix(rep(1/num.gauss,num.gauss),1,num.gauss)
  if (is.null(mu)) mu <- matrix(min(y)+diff(range(y))*seq(0,1,length=num.gauss+2)[-c(1,num.gauss+2)],1,num.gauss)
  if (is.null(sigma)) sigma <- matrix(rep(sd(y)/sqrt(num.gauss),num.gauss),1,num.gauss)

  llk <- logLik(y,p[1,],mu[1,],sigma[1,])
  while (flag) {
    count <- count+1
    # E step
    E <- matrix(0,n,num.gauss)
    p.s.k <- p[count,]/sigma[count,] # calculating p.k/sigma.k for all k
    for (i in 1:n) {
      exp.i <- exp(-((y[i]-mu[count,])^2)/(2*(sigma[count,]^2)))
      for (k in 1:num.gauss) {
        E[i,k] <- p.s.k[k]*exp.i[k]/sum(p.s.k*exp.i)
      }
    }
    rm(exp.i,p.s.k)

    # M step
    mu <- rbind(mu,rep(0,num.gauss))
    sigma <- rbind(sigma,rep(0,num.gauss))
    p <- rbind(p,rep(0,num.gauss))
    for (k in 1:num.gauss) {
      mu[count+1,k] <- sum(E[,k]*y)/sum(E[,k])
      sigma[count+1,k] <- sqrt(sum(E[,k]*(y-mu[count,k])^2)/sum(E[,k]))
      p[count+1,k] <- sum(E[,k])/(sum(E))
    }

    #m <- 0
    #temp <- sqrt(sum(((mu[count+1,]-mu[count,]))^2))
    #if (temp>m) m <- temp
    #temp <- sqrt(sum(((sigma[count+1,]-sigma[count,]))^2))
    #if (temp>m) m <- temp
    #temp <- sqrt(sum(((p[count+1,]-p[count,]))^2))
    #if (temp>m) m <- temp
    llk <- c(llk,logLik(y,p[count+1,],mu[count+1,],sigma[count+1,]))
    m <- abs(llk[count+1]-llk[count])


    if (m<=err) {
      flag <- F
    } else {
      if (count>=max.iter) {
        cat("Maximum iteration reached\n")
        flag <- F
      }
    }
  }
  return (list(mu=mu[count+1,],sigma=sigma[count+1,],p=p[count+1,],iterations=count,max.llk=max(llk)))
}


logLik <- function(y,p,mu,sigma) {
  v <- rep(0,length(y))
  for (i in 1:length(p)) {
    v <- v+p[i]*dnorm(y,mu[i],sigma[i])
  }
  llk <- log(v)
  min.obs <- min(llk[!is.infinite(llk)])
  llk[is.infinite(llk)] <- min.obs
  return(sum(llk))
}

# plots density function of mixture model
# finds intersection between each two gaussians
# which.intersect is the number of intersection to return
#                 (1 - between the first and second Gaussians from the left)
#                 (0 - intersection between the leftmost and the rightmost Gaussians)
# save.plot - path in which to save the density plot. If NULL => no save

calcIntersection <- function (y,obj=NULL,x.step=0.01,mu,sigma,p,draw.point=T,which.intersect=1,x.lab="log(v.max+2)",
                               save.plot=NULL,do.plot=T) {
  if (!is.null(obj)) {
    mu <- obj$mu
    sigma <- obj$sigma
    p <- obj$p
  }
  x <- seq(min(y)-0.25*diff(range(y)),max(y)+0.25*diff(range(y)),by=x.step)
  d <- matrix(0,length(x),length(mu))
  ncol.d <- ncol(d)
  if (do.plot) {
    if (!is.null(save.plot)) png(filename=save.plot,width=1000,height=1000)
    plot(density(y),main="",xlab=x.lab,ylab="density")
  }
  for (i in 1:ncol.d) {
    d[,i] <- dnorm(x,mean=mu[i],sd=sigma[i])*p[i]
    if (do.plot) lines(d[,i]~x,col="blue")
  }
  d.sum <- apply(d,1,sum)
  if (do.plot) lines(d.sum~x,col="red")
  if (length(mu)==1) return(numeric())
  w <- rep(0,ncol.d-1)
  x.intersect <- rep(0,ncol.d-1)
  y.intersect <- rep(0,ncol.d-1)
  for (i in 1:(ncol.d-1)) {
    a <- sigma[i]^2-sigma[i+1]^2
    b <- 2*(mu[i]*sigma[i+1]^2-mu[i+1]*sigma[i]^2)
    c <- 2*(sigma[i]^2)*(sigma[i+1]^2)*log(p[i]*sigma[i+1]/(sigma[i]*p[i+1]))+
      ((mu[i+1]^2)*(sigma[i]^2))-((mu[i]^2)*(sigma[i+1]^2))

    x.temp <- (-b+c(-1,1)*sqrt(b^2-4*a*c))/(2*a)
    x.temp <- x.temp[(x.temp>=mu[i])&(x.temp<=mu[i+1])]
    x.intersect[i] <- x.temp
    y.intersect <- dnorm(x.intersect[i],mu[i],sigma[i])*p[i]
    if ((draw.point)&&(do.plot)) points(x.intersect[i],y.intersect,pch=3,lwd=2)
  }
  # intersect between the first and last
  if (ncol.d>2) {
    a <- sigma[1]^2-sigma[ncol.d]^2
    b <- 2*(mu[1]*sigma[ncol.d]^2-mu[ncol.d]*sigma[1]^2)
    c <- 2*(sigma[1]^2)*(sigma[ncol.d]^2)*log(p[1]*sigma[ncol.d]/(sigma[1]*p[2]))+
      ((mu[ncol.d]^2)*(sigma[1]^2))-((mu[1]^2)*(sigma[ncol.d]^2))

    x.temp <- (-b+c(-1,1)*sqrt(b^2-4*a*c))/(2*a)
    x.temp <- x.temp[(x.temp>=mu[1])&(x.temp<=mu[ncol.d])]
    x.intersect[ncol.d] <- x.temp
    y.intersect <- dnorm(x.intersect[ncol.d],mu[ncol.d],sigma[ncol.d])*p[ncol.d]
    if ((draw.point)&&(do.plot)) points(x.intersect[ncol.d],y.intersect,pch=3,lwd=2)
  }
  if ((!is.null(save.plot))&&(do.plot)) dev.off()
  which.intersect <- pmin(ncol.d-1,which.intersect)
  if (which.intersect==0) return(x.intersect[length(x.intersect)])
  return(x.intersect[which.intersect])
}

