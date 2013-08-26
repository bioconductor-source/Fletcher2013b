#######################################################
###Customized heatmap function for large  sym. matrix
#######################################################
dot4heatmap <- function (x,px=NULL,xdendro = NULL,distfun = dist, 
                         hclustfun = hclust,reorderfun = function(d,w) reorder(d, w),add.expr, na.rm = TRUE,margins = c(5, 5), 
                         useRaster=FALSE,  ColSideColors, RowSideColors,LegSideColors,legwidth=0.25,sidewidth=0.15,
                         matwidth=3.5,dendwidth=0.6, breaks, cexRowCol = 0.2 + 1/log10(nrc), labRowCol = NULL, cexLab = 1.0, 
                         cexLeg=0.5,main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE,verbose = getOption("verbose"),
                         plotType = c("standard","points","userjob"),plotLayers = FALSE,pointSz = 1.0,plotScale = FALSE,
                         colorVec = NULL,zlim = c(-1,1),getOrd=FALSE,bgplot="white",qqCor=NULL,interCor=TRUE,interPlot=FALSE,
                         toCut=NULL,comment=NULL,...)
  {
  if( length(plotType) > 1 ) plotType <- "standard"
  if(plotType!="userjob"){
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
      stop("'x' must be a numeric matrix")
    nrc <- di[1L]
    if (nrc <= 1)
      stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2L)
      stop("'margins' must be a numeric vector of length 2")
    doDend <- !identical(xdendro, NA)      
    if (is.null(xdendro))
      xdendro <- rowMeans(x, na.rm = na.rm)
    if (doDend) {
      if (inherits(xdendro, "dendrogram"))
        dend <- xdendro
      else {
        hcr <- hclustfun(distfun(x))
        dend <- as.dendrogram(hcr)
        if (!is.logical(xdendro) || xdendro)
          dend <- reorderfun(dend, xdendro)
      }
      if (nrc != length(rcIndex <- order.dendrogram(dend)))
        stop("dendrogram ordering gave index of wrong length")
    }
    else rcIndex <- 1L:nrc
    x <- x[rcIndex, rcIndex]   
    colsIndex=rcIndex
    labRowCol <- if (is.null(labRowCol))
      if (is.null(rownames(x)))
        (1L:nrc)[rcIndex]
    else rownames(x)
    else labRowCol[rcIndex]     
  } else {
    if (is.null(labRowCol)){
      labRowCol<-1:length(unique(x[,1]))
    } 
    nrc<-length(unique(x[,1]))
    dend <- xdendro
    doDend=T
    rcIndex <- 1L:nrc	
    if (nrc != length(colsIndex<-order.dendrogram(dend)))stop("dendrogram ordering gave index of wrong length")
  }
  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(if (doDend) 1 else 0.05, 4)
  lhei <- c((if (doDend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 4)
  if (!missing(ColSideColors)) {
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1L], 0.2, lhei[2L])
  }
  if (!missing(RowSideColors)) {
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
    lwid <- c(lwid[1L], 0.2, lwid[2L])
  }
  lmat[is.na(lmat)] <- 0
  if (verbose) {
    cat("layout: widths = ", lwid, ", heights = ", lhei,"; lmat=\n")
    print(lmat)
  }
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  lwid =c(dendwidth,rep(sidewidth,4),matwidth)
  lhei =c(dendwidth,rep(sidewidth,4),matwidth)
  lmat<-rbind(c(0,0,0,0,0,11),c(12,0,0,0,0,8),c(12,0,0,0,0,7),c(13,0,0,0,0,6),c(13,0,0,0,0,5),c(10,4,3,2,1,9))
  layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
  if (!missing(RowSideColors)) {
    slabs=c("","","")
    if(!is.null(names(RowSideColors)))slabs=names(RowSideColors)
    par(mar = c(margins[1L], 0, 0.4, 0.2))
    image(rbind(1L:nrc), col = RowSideColors[[1]][colsIndex], axes = FALSE, useRaster=useRaster)
    mtext(side=1,slabs[1], las=2,cex=cexLeg)
    image(rbind(1L:nrc), col = RowSideColors[[2]][colsIndex], axes = FALSE, useRaster=useRaster)
    mtext(side=1,slabs[2], las=3,cex=cexLeg)
    image(rbind(1L:nrc), col = RowSideColors[[3]][colsIndex], axes = FALSE, useRaster=useRaster)
    mtext(side=1,slabs[3], las=3,cex=cexLeg)
    image(rbind(1L:nrc), col = RowSideColors[[4]][colsIndex], axes = FALSE, useRaster=useRaster)
    mtext(side=1,slabs[4], las=3,cex=cexLeg)       
  }
  if (!missing(ColSideColors)) {
    slabs=c("","","")
    if(!is.null(names(ColSideColors)))slabs=names(ColSideColors)
    par(mar = c(0.2, 0.4, 0, margins[2L]))
    image(cbind(1L:nrc), col = ColSideColors[[1]][colsIndex], axes = FALSE, useRaster=useRaster)
    mtext(side=2,slabs[1], las=1,cex=cexLeg)
    image(cbind(1L:nrc), col = ColSideColors[[2]][colsIndex], axes = FALSE, useRaster=useRaster)
    mtext(side=2,slabs[2], las=1,cex=cexLeg)
    image(cbind(1L:nrc), col = ColSideColors[[3]][colsIndex], axes = FALSE, useRaster=useRaster)
    mtext(side=2,slabs[3], las=1,cex=cexLeg)
    image(cbind(1L:nrc), col = ColSideColors[[4]][colsIndex], axes = FALSE, useRaster=useRaster)
    mtext(side=2,slabs[4], las=1,cex=cexLeg)        
  }
  par(mar = c(margins[1L], 0.3, 0.3, margins[2L]))
  if(plotType=="standard"){
    if(is.null(colorVec)) {
      colorVec<-heat.colors(12) 
    }
    if(!missing(breaks)){
      if(min(x)<breaks[1])breaks<-c(min(x), breaks)
      if(max(x)>breaks[length(breaks)])breaks<-c(breaks,max(x))
      colorVec<-colorRampPalette(colors=colorVec)(length(breaks)-1)	     	
    } else {
      colorVec<-colorRampPalette(colorVec, space="rgb")(20)
    }
    image(1L:nrc, 1L:nrc, x, xlim = 0.5 + c(0, nrc), ylim = 0.5 + c(0, nrc),breaks=breaks, 
          zlim = zlim, axes = FALSE, xlab = "", ylab = "", col = colorVec, useRaster=useRaster)
    axis(1, 1L:nrc, labels = labRowCol, las = 2, line = -0.5, tick = 0, cex.axis = cexRowCol)
    axis(4, 1L:nrc, labels = labRowCol, las = 2, line = -0.5, tick = 0, cex.axis = cexRowCol)         
  } else if(plotType=="points"){
    len <-length(x[,1])
    x <- matrix(data = x, nrow = prod(dim(x)), ncol = 1)
    xy <- rep(1:len, each = len)
    temp <- rep(1:len, len)
    xy<-cbind(xy,temp)
    xy<-xy[x>0,]
    x<-x[x>0]
    x<-cbind(xy,x)
    xy<-NULL
    if(is.null(colorVec)){
      colorVec<-c('red','yellow','black')
    }   
    if(!missing(breaks)){
      bkcenter=(breaks[-length(breaks)]+breaks[-1])/2
      bkcenter=c(-Inf,bkcenter,+Inf)
      cols=colorRampPalette(colors=colorVec)(length(bkcenter)-1)   
      x.col=rep(NA,length(x[,3]))
      cuts=cut(x[,3][!is.na(x[,3])], breaks=bkcenter,include.lowest=TRUE)
      x.col[!is.na(x[,3])]=cols[as.integer(cuts)]
      na.col="white"
      x.col[is.na(x.col)]=colorRampPalette(colors=c(na.col,na.col))(1)		
      colorVec=x.col	     	
    } else {
      temp <-as.factor(x[,3])
      colorVec <- colorRampPalette(colorVec, space="rgb")(length(levels(temp)))
      colorVec <- colorVec[temp]
    }
    plot(x[,1],x[,2], axes = FALSE, xlab = "", ylab = "",type="p",
         pch=15,col=colorVec,xlim = 0.5 + c(0, nrc), ylim = 0.5 + c(0, nrc),xaxs='i',yaxs='i',cex=pointSz)
    if (!plotScale){
      axis(1, 1L:nrc, labels = labRowCol, las = 2, line = -0.5, tick = 0,cex.axis = cexRowCol)
      axis(4, 1L:nrc, labels = labRowCol, las = 2, line = -0.5, tick = 0,cex.axis = cexRowCol)           
    } else {
      int<-round((nrc/3),0)
      axis(1, c(1,int,int*2, nrc),tick = T,line = 0.35,tcl=-0.2, mgp=c(3, 0.5, 0), cex.axis = cexRowCol )
      axis(4, c(1,int,int*2, nrc), tick = T,line = 0.3, las = 2,tcl=-0.2 ,mgp=c(3, 0.5, 0), cex.axis = cexRowCol)
    } 		
  } else if(plotType=="userjob"){          
    if(is.null(colorVec)){
      colorVec<-c('yellowgreen','yellow','red')
    }
    if(!is.null(toCut)){
      x <- x[ x[ ,3] > toCut , ]
    }
    xposit<-abs( x[,1] - x[,2] )
    xcor<-(x[,3])
    if(is.null(qqCor)){                                                
      temp <-as.factor(xcor) 
    } else {                      
      if(interCor){
        rankx<-abs(rank(xcor)+rank(1/xposit))
      } else {
        rankx<-rank(xcor)
      }
      #
      grNumber<- 1/qqCor;
      qx=quantile(rankx, probs = seq(0,1,grNumber) )         
      rankx <- cut(rankx, qx, labels = FALSE, include.lowest = TRUE)          
      temp <-as.factor(rankx)       
    }       
    colorVec <- colorRampPalette(colorVec, space="rgb")(length(levels(temp)))
    p.color <- colorVec[temp]
    if(!plotLayers){
      plot(1, axes = FALSE, xlab = "", ylab = "",type="n",
           pch=15,col=p.color,xlim = 0.5 + c(0, nrc), ylim = 0.5 + c(0, nrc ),
           xaxs='i',yaxs='i',cex=pointSz)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], border = NA, col = bgplot)               
      points(x[ ,1:2], type="p",pch=15,col=p.color,cex=pointSz)                  
    } else {
      #plota apenas matriz
      plot(1, axes = FALSE, xlab = "", ylab = "",type="n",
           pch=15,col=p.color[temp],xlim = 0.5 + c(0, nrc), 
           ylim = 0.5 + c(0, nrc ),xaxs='i',yaxs='i',cex=pointSz)               
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], border = NA, col = bgplot)
      #Define ranking de plotagem. Pode interpolar (cor - posicao) ordenando pelo menor dx desses dois rankings. 
      #Resultado e a interpolacao desses ranks. Do contrario, ordem de plotagem segue somente camadas de cor.
      if(interPlot){
        rankx<-abs(rank(xcor)-rank(xposit))
      } else {
        rankx<-rank(xcor)
      }
      #ordenar por pelo dx 
      indx<-sort.list(rankx)
      #plota por camadas segundo dx
      points(x[ indx ,1:2], type="p",pch=15,col=p.color[indx],cex=pointSz)                                              
    }
    if (!plotScale){
      axis(1, 1L:nrc, labels = labRowCol, las = 2, line = -0.5, tick = 0,cex.axis = cexRowCol)
      axis(4, 1L:nrc, labels = labRowCol, las = 2, line = -0.5, tick = 0,cex.axis = cexRowCol)            
    } else {
      int<-round((nrc/3),0)
      axis(1, c(1,int,int*2, nrc),tick = T,line = 0.35,tcl=-0.2, mgp=c(3, 0.5, 0))
      axis(4, c(1,int,int*2, nrc), tick = T,line = 0.3, las = 2,tcl=-0.2 ,mgp=c(3, 0.5, 0))
    }       
  }
  if(!is.null(px) & plotType=="standard"){
    a<-match(px[,1],colnames(x))
    b<-match(px[,2],rownames(x))
    xx<-a
    yy<-b
    idx<-a<b   
    xx[idx]<-b[idx]
    yy[idx]<-a[idx]
    points(x=xx,y=yy, pch="*",cex=1.5,col="black")
  } 
  if (!is.null(xlab)){
    mtext(xlab, side = 1, line = margins[1L] - 1.25,cex=cexLab)
  }
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2L] - 1.25, cex=cexLab)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  par(mar = c(margins[1L]+0.6, 0, 0.4, 0.1))
  if (doDend)
    plot(dend, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  else frame()
  par(mar = c(0.1, 1.0, if (!is.null(main)) 1 else 0.2, margins[2L]))
  if (doDend)
    plot(dend, axes = FALSE, xaxs = "i", leaflab = "none")
  else if (!is.null(main))
    frame()
  if (!is.null(main)) {
    par(xpd = NA)
    title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  if (!missing(LegSideColors)) {
    legwidth=0.5-legwidth
    slabs=c("","")
    if(!is.null(names(LegSideColors)))slabs=names(LegSideColors)
    par(mai=c(legwidth, 0.03, 0, 0))
    nrc=length(LegSideColors[[1]])
    image(cbind(1L:nrc), col = LegSideColors[[1]], axes = FALSE)
    axlabels=TRUE
    at=NULL
    if(!is.null(names(LegSideColors[[1]]))){
      axlabels=names(LegSideColors[[1]])
      at=seq(0,1, length.out=length(axlabels))
    }
    axis(side=1,at=at,cex.axis=cexLeg-0.1, tcl=-0.2, mgp=c(0,0,0),lwd=0.4,labels=axlabels);box(lwd=0.4)
    mtext(side=3,slabs[1], las=1,cex=cexLeg, adj=0)
    par(mai=c(legwidth/2, 0.03, legwidth/2, 0))
    nrc=length(LegSideColors[[2]])
    image(cbind(1L:nrc), col = LegSideColors[[2]], axes = FALSE)
    axlabels=TRUE
    at=NULL
    if(!is.null(names(LegSideColors[[2]]))){
      axlabels=names(LegSideColors[[2]])
      at=seq(0,1, length.out=length(axlabels))
    }
    axis(side=1, at=at,cex.axis=cexLeg-0.1, tcl=-0.2, mgp=c(0,0,0),lwd=0.4,labels=axlabels);box(lwd=0.4)
    mtext(side=3,slabs[2], las=1,cex=cexLeg, adj=0)
  }
  invisible(list(rcIndex = rcIndex, xdendro = if (keep.dendro && doDend) dend, 
                 xdendro = if (keep.dendro && rcIndex) dend))
  if(!is.null(comment))mtext(text=comment,side=1, cex=0.6, line = -1, outer=TRUE, adj=0)       
}


