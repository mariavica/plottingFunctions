##########
get.chr <- function (x, sep=sep, pos=1) {
  
   strsplit(x,split=sep)[[1]][pos]
  
}





correctAB <- function (x) {
  
  chrs <- unique(x[,1])
  
  for (i in 1:length(chrs)) {
    sel <- which(x[,1]==chrs[i])
    
    for (j in 5:ncol(x)) {
      
      if (cor(x[sel,4],x[sel,j]) < 0 ) {
        x[sel,j]<- -c(x[sel,j])
      }
      
      
    }
    
    
  }
  
  return(x)
  
}




selectAB <- function (x.list, ref.at) {
  
  x2.list<-list()
  
  for (zz in 1:length(x.list)) {
    
    x<-x.list[[zz]]
    x2 <- x[,1:4]
    x2[,4]<-0
    
    chrs <- unique(x[,1])
    
    id.x <- paste(x[,1],x[,2],x[,3],sep=':')
    id.ref <- paste(ref.at[,1],ref.at[,2],ref.at[,3],sep=':')
    
    rownames(x)<-id.x
    rownames(x2)<-id.x
    rownames(ref.at)<-id.ref
    
    for (i in 1:length(chrs)) {
      
      cat(paste(chrs[i],'\n'))
      
      id.x.sel <- id.x[which(x[,1]==chrs[i])]
      id.ref.sel <- id.ref[which(ref.at[,1]==chrs[i])]
      
      
      sel <- intersect(id.x.sel, id.ref.sel)
      
      if (length(sel)<=0) {warning('No common IDs, skipping chromosome...'); next}
      
      # correlation
      
      cor.1 <- cor(x[sel,4],ref.at[sel,4])
      cor.2 <- cor(x[sel,5],ref.at[sel,4])
      
      cors <- c(cor.1,cor.2)
      # pick the maximum
      sel.cor <-which(abs(cors)==max(abs(cors)))
      
      x2[sel,4] <- x[sel,sel.cor+3]*sign(cors[sel.cor])
      

			cat(paste('  Correlation with AT content: ', round(cors[1],2),', ', round(cors[2],2),'\n' ,sep=''))
			cat(paste('  Selected column: EV', sel.cor,'\n',sep=''))
			cat(paste('  Invert sign: ', sign(cors[sel.cor])==(-1),'\n',sep=''))
			cat('\n')

		}

		x2.list[[zz]]<-x2

	}

	return(x2.list)

}











plotTwoTrends <- function(x1, x2, x = NA, 
                          col1 = 'red', col2 = 'blue', lwd = 3,
                          normalization = 'none', norm_independent = FALSE,
                          expand = 0.25, fix.ylim = NULL,
                          loess = FALSE, sp = 0.5, res.x = 1000, 
                          plot.sd = TRUE, error.plot = 0.95,
                          xlab = NA, ylab1 = NA, ylab2 = NA, omit.lab2 = FALSE,
                          main = NULL, pos.leg='topright', name.1=NULL, name.2=NULL) {


  if (ncol(x1)!=ncol(x2)) {
    stop('X1 and X2 have different number of columns')
  }
  
  error.plot <- qnorm(1-(1-error.plot)/2)
  
  
  if (is.na(x)) { x <- 1:ncol(x1) }
 
  y1 <- as.vector(unlist(x1))
  y1.x <- rep(x,each=nrow(x1))
  
  if (loess) {
    if (length(y1)>10^7) {
      warning('X1 too large, taking 10^7 random numbers')
      sel <- sample(1:length(y1),10^7)
      y1 <- y1[sel]
      y1.x <- y1.x[sel]
    }
    
    lo <- loess(y1~y1.x,span=sp)
 

    xl <- seq(from = min(y1.x), to= max(y1.x), length.out = res.x)
    yl <- predict(lo,xl)
    
    
    if (is.null(fix.ylim)) {
      fix.ylim <-c(min (yl) -(max(yl)-min(yl))*expand , max(yl) + (max(yl)-min(yl))*expand )
    } 
  

    plot(xl, yl , col=col1, lwd=lwd, type='l', main=main, axes=FALSE, xaxt='n',ylab='', xlab='',
         ylim=fix.ylim)
    
    if (is.na(ylab1)) {ylab1 <- colnames(x1)}
  
    axis(side = 1 , at=1:ncol(x1), labels=ylab1)
    axis(side = 2, col=col1)
    
    
    if (plot.sd) {
      
      lox<- predict(loess(y1~y1.x,span=sp), newdata = xl ,se=TRUE)
      upper<- lox$fit + error.plot*lox$se.fit
      lower<- lox$fit - error.plot*lox$se.fit

      polygon( x = c(xl, rev(xl)), y=c(upper,rev(lower)),
               col=rgb(red=1,blue=0,green=0,alpha=0.25), border=NA)

    }
    
        
  } else {
    
    #regular plot
    
    yl <- tapply(y1,factor(y1.x), FUN=mean)
    xl <- as.numeric(unique(factor(y1.x)))

    if (is.null(fix.ylim)) {
      fix.ylim <-c(min (yl) -(max(yl)-min(yl))*expand , max(yl) + (max(yl)-min(yl))*expand )
    } 
    
    plot(xl, yl , col=col1, lwd=lwd, type='l', main=main, axes=FALSE, xaxt='n',ylab='', xlab='',
         ylim=fix.ylim)
    
    if (is.na(ylab1)) {ylab1 <- colnames(x1)}
    
    axis(side = 1 , at=1:ncol(x1), labels=ylab1)
    axis(side = 2, col=col1)
    
    
    if (plot.sd) {
      
      yl.sd <- tapply(y1,y1.x, FUN=sd)
      upper<- yl + error.plot*yl.sd
      lower<- yl - error.plot*yl.sd
      
      polygon( x = c(xl, rev(xl)), y=c(upper,rev(lower)),
               col=rgb(red=1,blue=0,green=0,alpha=0.25), border=NA)
      
    }
    
    
  }
 

 
  
  ## second plot
  
  y2 <- as.vector(unlist(x2))
  y2.x <- rep(x,each=nrow(x2))
  
  if (loess) {
    if (length(y2)>10^7) {
      warning('X2 too large, taking 10^7 random numbers')
      sel <- sample(1:length(y2),10^7)
      y2 <- y2[sel]
      y2.x <- y2.x[sel]
    }
    
    lo <- loess(y2~y2.x,span=sp)
    
    xl <- seq(from = min(y2.x), to= max(y2.x), length.out = res.x)
    yl <- predict(lo,xl)
    
    
    if (is.null(fix.ylim)) {
      fix.ylim <-c(min (yl) -(max(yl)-min(yl))*expand , max(yl) + (max(yl)-min(yl))*expand )
    } 
    
    par(new = TRUE)

    plot(xl, yl , col=col2, lwd=lwd, type='l', main=NULL, axes=FALSE, xaxt='n',ylab='', xlab='',
         ylim=fix.ylim)
    
    if (is.na(ylab2)) {ylab2 <- colnames(x2)}
    
    axis(side = 4, col=col2)
    
    if (plot.sd) {
      
      lox<- predict(loess(y2~y2.x,span=sp), newdata = xl ,se=TRUE)
      upper<- lox$fit + error.plot*lox$se.fit
      lower<- lox$fit - error.plot*lox$se.fit
      
      polygon( x = c(xl, rev(xl)), y=c(upper,rev(lower)),
               col=rgb(red=0,blue=1,green=0,alpha=0.25), border=NA)
      
    }
    

  } else {
    
    # regular plot
    yl <- tapply(y2,factor(y2.x), FUN=mean)
    xl <- as.numeric(unique(factor(y2.x)))
    
    if (is.null(fix.ylim)) {
      fix.ylim <-c(min (yl) -(max(yl)-min(yl))*expand , max(yl) + (max(yl)-min(yl))*expand )
    } 
    
    par(new = TRUE)
    
    plot(xl, yl , col=col2, lwd=lwd, type='l', main=NULL, axes=FALSE, xaxt='n',ylab='', xlab='',
         ylim=fix.ylim)
    
    if (is.na(ylab2)) {ylab2 <- colnames(x2)}
    
    axis(side = 4, col=col2)
    
    if (plot.sd) {
      
      yl.sd <- tapply(y2, y2.x, FUN=sd)
      upper<- yl + error.plot*yl.sd
      lower<- yl - error.plot*yl.sd
      
      polygon( x = c(xl, rev(xl)), y=c(upper,rev(lower)),
               col=rgb(red=0,blue=1,green=0,alpha=0.25), border=NA)
      
    }
    
    
    
    
  }
  
  
  ### Legend
  
  if(is.null(name.1)) deparse(substitute(x1))
  if(is.null(name.2)) deparse(substitute(x2))
  

  legend(pos.leg, c(name.1,name.2), col=c(col1,col2), lwd=2, bty = "n")
  
  
  
  
}

  
  
 # a <- matrix(rnorm(1000), ncol=10, nrow = 100)
#  b <- matrix(rnorm(1000), ncol=10, nrow = 100)
  
 # plotTwoTrends(x1 = a, x2 = b , loess = TRUE , expand = 0.75)
  
  plotTwoTrends(x1 = a, x2 = b , loess = TRUE , expand = 0.75, main = 'Juju', name.1 = 'a', name.2 = 'b')
  











  library(h5)
  testmat <- matrix(rpois(1000000,lambda=50), ncol = 1000, nrow=1000)
  # Create HDF5 File
  file <- h5file("test.h5")
  # Save matrix to file in group 'testgroup' and datasetname 'testmat'
  file["testgroup", "testmat"] <- testmat
  # Close file
  h5close(file)
  
  
  
  
  
  
  
  
  
  
  