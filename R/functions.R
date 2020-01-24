##########

#####

#ref <- data.frame(chr=c('chr1'),start=c(300,600), end=c(400,700))

#map <- data.frame(chr=rep('chr1',10),start=c(100,200,300,400,500,600,700,800,900,1000),
#                  end=c(100,200,300,400,500,600,700,800,900,1000)+99,
#                  V1=rnorm(10),
#                  V2=rnorm(10,mean = 2))


#map <- read.table('PCM.index.bed', header = FALSE, sep='\t')
#ref <- read.table('fc3e8b36a_dec8ab5e3.borders.bed', header = FALSE, sep='\t')



map.plot <- function ( ref, map , grouping, expand=5, show.orig=FALSE ) {
  
  # check columns
  ref[,1]<-as.character(ref[,1])
  map[,1]<-as.character(map[,1])
  
  ref[,2]<-as.numeric(ref[,2])
  map[,2]<-as.numeric(map[,2])
  
  ref[,3]<-as.numeric(ref[,3])
  map[,3]<-as.numeric(map[,3])
  
  colnames(map)[c(1,2)]<-c('chr','start')
  
  
  #sort map
  map <- map[with(map, order(chr, start)),]
  
  select.all <- matrix(0, ncol=(ncol(map)-3), nrow=expand*2+1)

  numbers <- matrix(0, ncol=(ncol(map)-3), nrow=expand*2+1)
  
  
  for (i in 1:nrow(ref)) {

    # find 
    sel <- which(map[,1]==ref[i,1])
    
    sel.more <- which(map[sel,2] <= ref[i,2] &  map[sel,2] >= ref[i,2] )  # identify the right thing
    
    if (length(sel.more)==1) {
    
      sel.sel <- c( (sel.more-expand) : (sel.more+expand) )
      sel.sel[which(sel.sel<0)]<-0
      
      tocp <-map[ sel[ sel.sel   ] , 4:ncol(map) ]
      
      numbers <- numbers + as.numeric(!is.na(tocp))
      
      for (ii in 1:nrow(tocp)) {
        for (jj in 1:ncol(tocp)) {
          if (is.na(tocp[ii,jj])==TRUE) tocp[ii,jj]<-0
        }
      }
      
      select.all <- select.all + tocp
      
    }

  }
  
  select.all <- select.all/numbers
  
  ### smooth 
  
  plot (  spline(1:(expand*2+1), select.all[,1], n=1000 ), type='l', ylim=c(min(select.all),max(select.all)  ),
          #main=paste(coordinates[1,1], coordinates[1,2], coordinates[nrow(coordinates),3], sep=':'),
          ylab='signal',xlab='position',xaxt='n')

  step <- map[1,3]-map[1,2]
  
  axis(1, at=1:(expand*2+1), c(  rev( -(1:expand)*step ), 0,  (1:expand)*step ) , las=2 )

  for (j in 2:ncol(select.all)) {
    
    lines(spline(1:(expand*2+1), select.all[,j], n=1000 ), col=j)
    
  }
  
  
  if (show.orig) {
    
    for (j in 1:ncol(select.all)) {
      lines(1:(expand*2+1), select.all[,j], col=j, lty=2)
    }

  }

  
}


#map.plot(ref=ref,map=map, expand=5, show.orig = TRUE)







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



GOheatmap <- function(genes, cluster=NA, plot=TRUE, organism="mouse", minGSSize = NA,
                      maxGSSize = NA, colnames = NA, cutoff.top=NA, cutoff.dw=NA,
                      simplify.cutoff = NA, labels_col=unique(cluster), fdr.cutoff = 1) {
  
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package \"clusterProfiler\" needed for this function to work. Please install it.",
         call. = FALSE)
  } 
  
  #lab <- 
  
  require(clusterProfiler)
  
  #organism="human"
  
  if (organism=="human") { require(org.Hs.eg.db); organism<-org.Hs.eg.db }
  if (organism=="mouse") { require(org.Mm.eg.db); organism<-org.Mm.eg.db }
  
  #require(org.Hs.eg.db); organism<-org.Hs.eg.db
  
  if (is.list(genes)) {
 
  } else {
    
  }
  
  #get all possible GO terms
  results <- enrichGO(gene         = as.character(genes),
                      OrgDb         = organism,
                      keyType       = 'SYMBOL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize,
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1)
  
  if (!is.na(simplify.cutoff)) {
    results <- simplify(results,cutoff=simplify.cutoff)
  }
  
  mat.pval <- matrix(1, ncol=length(unique(cluster)), nrow= nrow(results) )
  rownames(mat.pval) <- results$ID
  description <- results$Description

  mat.padjust <- mat.pval
  mat.generatio <- mat.pval
  mat.qvalue <- mat.pval
  
  #analyze cluster by cluster
  for (i in 1:length(unique(cluster))) {
    
    results <- enrichGO(gene         = as.character(genes)[which(cluster==i)],
                        OrgDb         = organism,
                        keyType       = 'SYMBOL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        minGSSize = minGSSize,
                        maxGSSize = maxGSSize,
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1)
    sel <- which((results$ID %in% rownames(mat.pval)) == TRUE)
    
    mat.pval[results$ID[sel],i]<-results$pvalue[sel]
    mat.padjust[results$ID[sel],i]<-results$p.adjust[sel]
    mat.generatio[results$ID[sel],i]<-results$GeneRatio[sel]
    mat.qvalue[results$ID[sel],i]<-results$qvalue[sel]
    print(i)
    
  }  
    
    
  require(pheatmap)
  #print if needed
  if (plot) {
    
    sel <- which(apply(mat.padjust,1,min)<fdr.cutoff)
    
    toplot <- -log10(mat.padjust[sel,])
    
    print(head(toplot))
    
    toplot2<-toplot
    
    if (!is.na(cutoff.top)) {
      for (i in 1:nrow(toplot)) {
        for (j in 1:ncol(toplot)) {
          if (toplot[i,j]>cutoff.top) toplot[i,j]<-cutoff.top
        }
      }
    }
    
    print(head(toplot))

    if (!is.na(cutoff.dw)) {
      for (i in 1:nrow(toplot)) {
        for (j in 1:ncol(toplot)) {
          if (toplot[i,j]<cutoff.dw) toplot[i,j]<-cutoff.dw
        }
      }
    }
    
    
        
    pheatmap(toplot, clustering_distance_rows='correlation', 
             cluster_cols = FALSE, labels_row = description[sel],
             labels_col=labels_col, color=colorRampPalette(c("white", "blue"))( 20))
    
    
  } else {

   

  }
  
  return(list(   mat.pval, mat.padjust, mat.generatio, mat.qvalue, description[sel]
  ))
  
  
}




plotStackedGenes <- function (x, groups) {
  
  means <- aggregate(t(x), by=list(groups), mean)
  
  shap <- barplot(t(as.matrix(means[,-1])),beside=TRUE, ylim=c(-2,2))
  
  tt <- table(groups)
  
  for (i in 1:nrow(x)) {
    
    for (j in 1:nrow(means)) {
      
      points(rep(shap[i,j], tt[j]), x[i,(1:tt[j])+sum(tt[(0):(j-1)])], pch=19)

    }
  
    
  }
  
  
  
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



writeGsea <- function (x, cath, filename) {

  # gmt
  sink(file=paste(filename,"_head.gmt",sep=''))
  cat('#1.2\n')
  cat(paste(nrow(x),"\t",ncol(x),'\n',sep=''))
  sink()

  x.p<-data.frame(NAME=rownames(x),Description=rownames(x),x)
  write.table(x.p,paste(filename,"_body.gmt",sep=''),quote=FALSE,row.names=FALSE)
  
  system(paste('cat ',filename,'_head.gmt ',filename,"_body.gmt > ",filename,'.gmt',sep=''))
  system(paste('rm ',filename,'_head.gmt',sep=''))
  system(paste('rm ',filename,'_body.gmt',sep=''))
  
  #* cls
  
  if (length(cath)==1) {
    
    ### cls is the parameter you've entered (numerical category)
    sink(paste(filename,".cls",sep=''))
    cat("#numeric\n")
    cat(paste('#',cath,'\n',sep=''))
    cat(unlist(x[cath,]))
    sink()

  } else {
 
    ## cls is a factor category

  }

}


#jaja <- matrix(1:10,ncol=2,nrow=5)
#rownames(jaja)<-letters[1:5]
#colnames(jaja)<-LETTERS[1:2]

#write.gsea(jaja,cath='a',filename = "juju")


  
  
  
  
  
  
  