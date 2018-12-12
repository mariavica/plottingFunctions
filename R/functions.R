
plotTwoTrends <- function(x1, x2, x = NA, 
                          col1 = 'blue', col2 = 'red', lwd = 3,
                          normalization = 'none', norm_independent = FALSE,
                          expand = 0.25, fix.ylim = NULL,
                          loess = FALSE, sp = 0.5, res.x = 1000, 
                          plot.sd = TRUE, 
                          xlab = NA, ylab1 = NA, ylab2 = NA, omit.lab2 = FALSE ) {


  if (ncol(x1)!=ncol(x2)) {
    stop('X1 and X2 have different number of columns')
  }
  
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
 
    sds<-loess.sd(y1~y1.x, nsigma = 1, span=sp)

    xl <- seq(from = min(y1.x), to= max(y1.x), length.out = res.x)
    yl <- predict(lo,xl)
    
    
    if (is.null(fix.ylim)) {
      ylim.toplot <-c(min (yl) -(max(yl)-min(yl))*expand , max(yl) + (max(yl)-min(yl))*expand )
    } 
    

    plot(xl, yl , col=col1, lwd=lwd, type='l', main=NULL, xaxt='n',ylab='', xlab='',
         ylim=fix.ylim)
    
    if (is.na(ylab1)) {ylab1 <- colnames(x1)}
  
    axis(side = 1 , at=1:ncol(x1), labels=ylab1)
 
    if (plot.sd) {
      
      lox<- predict(loess(y1~y1.x,span=sp), newdata = xl ,se=TRUE)
      upper<- lox$fit + 2*lox$se.fit
      lower<- lox$fit - 2*lox$se.fit

      polygon( x = c(xl, rev(xl)), y=c(upper,rev(lower)),
               col=rgb(red=1,blue=0,green=0,alpha=0.25), border=NA)

    }
    
        
  } else {
    
    
    plot()
    
    
  }
 

 
  
  
  
  xl <- seq(min(hic.x),max(hic.x), (max(hic.x) - min(hic.x))/1000)
  yl <- predict(lo,xl)
  plot(xl, yl , col='red', lwd=3, type='l', main=NULL, xaxt='n',ylab='', xlab='',
       ylim=c(min (yl) - (max(yl)-min(yl))*expand , max(yl) + (max(yl)-min(yl))*expand )
       #ylim=c(-0.06, 0.04)
  )
  
  
  time.2<-seq(from=1,to=6,length.out = 1000)
  
  
  lox<- predict(loess(hic.y~hic.x,span=sp), newdata = time.2 ,se=TRUE)
  upper<- lox$fit + 2*lox$se.fit
  lower<- lox$fit - 2*lox$se.fit
  
  
  
  polygon( x = c(time.2, rev(time.2)), y=c(upper,rev(lower)), col=rgb(red=1,blue=0,green=0,alpha=0.25), border=NA)
  
  sel<-to.cluster[which(all$cluster==j),1:12]  #check
  hic.cor.y<-apply( sel, 2, mean)
  hic.loess.y<-lox$fit
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 
 
  
  
  
  time<-1:6
  
  hic.y<-as.vector(unlist(to.cluster[which(all$cluster==j),1:12]))    ### check
  hic.x<-rep(rep(time,each=2),each=length(which(all$cluster==j)))          ### check
  
  lo <- loess(hic.y~hic.x,span=sp)
  
  
  sds<-loess.sd(hic.y~hic.x, nsigma = 1, span=sp)
  
  #plot(x,y)
  xl <- seq(min(hic.x),max(hic.x), (max(hic.x) - min(hic.x))/1000)
  yl <- predict(lo,xl)
  plot(xl, yl , col='red', lwd=3, type='l', main=NULL, xaxt='n',ylab='', xlab='',
       ylim=c(min (yl) - (max(yl)-min(yl))*expand , max(yl) + (max(yl)-min(yl))*expand )
       #ylim=c(-0.06, 0.04)
  )
  
  axis(side = 1 , at=1:6, labels=c("0h","1h","3h","6h","18h","120h"))
  
  time.2<-seq(from=1,to=6,length.out = 1000)
  
  
  lox<- predict(loess(hic.y~hic.x,span=sp), newdata = time.2 ,se=TRUE)
  upper<- lox$fit + 2*lox$se.fit
  lower<- lox$fit - 2*lox$se.fit
  
  
  
  polygon( x = c(time.2, rev(time.2)), y=c(upper,rev(lower)), col=rgb(red=1,blue=0,green=0,alpha=0.25), border=NA)
  
  sel<-to.cluster[which(all$cluster==j),1:12]  #check
  hic.cor.y<-apply( sel, 2, mean)
  hic.loess.y<-lox$fit
  
  
  ## rna
  
  rna.y<-as.vector(unlist(to.cluster[which(all$cluster==j),13:24]))  #check
  rna.x<-rep(rep(time,2),each=length(which(all$cluster==j)))       #check
  
  
  par(new = TRUE)
  
  lo <- loess(rna.y~rna.x,span=sp)
  #plot(x,y)
  xl <- seq(min(rna.x),max(rna.x), (max(rna.x) - min(rna.x))/1000)
  yl <- predict(lo,xl)
  plot(xl, yl, col='blue', lwd=3, type='l',axes=FALSE,xaxt='n',ylab='',
       ylim=c(min (yl) - (max(yl)-min(yl))*expand , max(yl) + (max(yl)-min(yl))*expand )
       #ylim=c(-0.25,1)
  )
  if (i %in% good) { 
    
    
    lox<- predict(lo, newdata = time.2 ,se=TRUE)
    upper<- lox$fit + 2*lox$se.fit
    lower<- lox$fit - 2*lox$se.fit
    
    time.2<-seq(from=1,to=11,length.out = 1000)
    
    
    polygon( x = c(time.2, rev(time.2)), y=c(upper,rev(lower)), col=rgb(red=0,blue=1,green=0,alpha=0.25), border=NA)
    
  }
  sel<-to.cluster[which(all$cluster==j),13:24]   ### check
  rna.cor.y<-apply( sel, 2, mean)
  rna.loess.y<-lox$fit
  
  
  cor.data[j]<-cor(hic.cor.y,rna.cor.y, method = 'spearman')
  p.data[j]<-cor.test(hic.cor.y,rna.cor.y, method = 'spearman')$p.val
  cor.loess[j]<-cor(hic.loess.y,rna.loess.y, method = 'spearman')
  p.loess[j]<-cor.test(hic.loess.y,rna.loess.y, method = 'spearman')$p.val
  
  
  title(main=paste( "Cluster ",j," (cor=",round(cor.data[j],3),', p=',round(p.data[j],3),')',sep="")  )
  legend("bottomright",c("WT","R35A"),col=c("red","blue"),lwd=3)     ##### check
  
  axis(side = 4)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}
