##################################################################################
### Supplementary functions only used to generate  main figs. for Fletcher et al.
##################################################################################
#function to plot aggrement among signatures
plotcon<-function(L1,tlt="MRA agreement (Cohort I)",
                  sig=1e-3, top=NULL, glabs=NULL, cols=NULL, vennsz=0.8,
                  cexplab=0.75, cexpaxis=1.0, cexp=0.3, cexdiag=1.35, cexvenn1=1.0,
                  cexvenn2=0.7, cex.axis=1.0, labp="Regulons (enrichment rank)",
                  lfaxispos=1,btaxispos=1,lwd=1.0, tcl=-0.5){
  llabs<-names(L1)
  rlabs<-rownames(L1)
  if(is.null(glabs))glabs<-names(L1)
  ngroups=ncol(L1)
  if(is.null(cols)){
    cols<-rainbow(ngroups, alpha = 1)
  }
  maxg<-0
  for(i in 1:ncol(L1)){
    pt<-sum(!is.na(L1[,i]))
    if(pt> maxg)maxg=pt
  }
  narmv<-function(x,y,labs){
    na.map<-is.na(x)+is.na(y)
    idx<-na.map==0
    xy<-cbind(x[idx],y[idx])
    rownames(xy)<-labs[idx]
    colnames(xy)<-c("x","y")
    xy  	
  }
  getRanks<-function(x,y,labs){
    xy<-narmv(x,y,labs)
    xy[,1]<-rank(xy[,1], na.last=NA )
    xy[,2]<-rank(xy[,2], na.last=NA )		
    xy
  }
  vc<-matrix(NA, nrow=(ngroups*(ngroups-1))/2, ncol=2)
  n=1
  for(i in 1:ngroups){
    for(j in 1:ngroups){
      if(i<j){
        vc[n,]<-c(i,j)
        n=n+1				
      }
    }
  }
  btaxispos=n-btaxispos
  sizes<-list()
  bgsz=0
  for(i in 1:nrow(vc)){
    if(is.null(top)){
      x=L1[,vc[i,1]]
      y=L1[,vc[i,2]]
      xy=narmv(x,y, rlabs)
      px=p.adjust(xy[,1])<sig
      py=p.adjust(xy[,2])<sig			
    } else {
      x=L1[,vc[i,1]]
      y=L1[,vc[i,2]]
      xy=getRanks(x,y, rlabs)
      px=xy[,1]<=top
      py=xy[,2]<=top			
    }
    labs<-rownames(xy)	
    szi<-length(intersect(labs[px],labs[py]))
    szu<-length(union(labs[px],labs[py]))
    sza<-sum(px)
    szb<-sum(py)
    if(sza>bgsz)bgsz=sza
    if(szb>bgsz)bgsz=szb
    invert=ifelse(sza>=szb,1,2)	
    sizes[[i]]<-c(sza, szb, szi,szu)
  }
  j=1
  panel.venn<-function(x, y,...){		
    sz<-sizes[[j]]
    relsz=(sz[1:2]/bgsz)/10
    invert=ifelse(sz[1]>=sz[2],1,2)
    pp<-draw2venn(sz[1], sz[2], sz[3], inverted=c(FALSE,TRUE)[invert], max.circle.size=sqrt(max(relsz))*vennsz)
    xlim <- ylim <- c(0, 4)
    plot.window(xlim, ylim)
    ord=c(1:2)
    posx=c(0.4,3.6)
    if(invert==2)ord=rev(ord)
    for(i in ord){
      xx=pp[[i]]$x*4
      yy=pp[[i]]$y*4
      cc=adjustcolor(cols,alpha.f=0.5)
      polygon(x=xx,y=yy,col=cc[vc[j,i]],lty="dotted")
      text(x=posx[i],y=3.2,paste("n=", sz[i],sep=""),cex=cexvenn1)
    }
    jc=round(sz[3]/sz[4],2)
    mtext(side=1,paste("JC=",jc,sep=""), line=-1.5, cex=cexvenn2)
    xy=getRanks(x,y,rlabs)
    r <-cor.test(xy[,1],xy[,2], method="spearman", exact=FALSE)
    mtext(side=1,paste("R=",round(r$estimate,2),sep=""), line=-0.5, cex=cexvenn2)
    #pp<-format(r$p.value, scientific=TRUE, digits=2)
    #mtext(side=1,paste("P=",pp,sep=""), line=0.5, cex=0.6)			
    assign("j",value=get("j")+1, inherits=TRUE)
  }
  k=1
  panel.points<-function(x, y,...){
    xy=getRanks(x,y,rlabs)
    x=xy[,1]
    y=xy[,2]
    xlim <- c(0, maxg)
    ylim <- c(0, maxg)
    plot.window(xlim, ylim)
    points(x,y, pch=15, cex=cexp,lwd=0.5)
    #smoothScatter(x,y, bandwidth=12, add=TRUE)
    box(col = "lightgray")
    if(k%in%btaxispos){
      axis(1, cex.axis=cex.axis, lwd=lwd, tcl=tcl)
      mtext(labp,side=1,col="black",line=2.5,cex=cexplab)
    }
    if(k%in%lfaxispos){
      axis(2,cex.axis=cex.axis,lwd=lwd, tcl=tcl)
      mtext(labp,side=2,col="black",line=2.5,cex=cexplab)
    }
    assign("k",value=get("k")+1, inherits=TRUE)			
  }
  l=1
  textpanel<-function(...){
    plot.window(c(0,1), c(0,1))
    text(x=0.5,y=0.5,glabs[l],cex=cexdiag,col=cols[l])
    assign("l",value=get("l")+1, inherits=TRUE)
  }
  corrgram(L1, row1attop=TRUE, main=tlt, 
           gap=0, labels=glabs, lower.panel=panel.points, 
           upper.panel=panel.venn, text.panel=textpanel)
}
##-----------------------------------------------------------------------------
#simple function to cross labels in aggrement plots
#regulons in the agreement matrix are ether assigned to random values or NA, 
#and then updated with the observed values in each experiment.
getmatch<-function(lts,glabs, set.seed.na=FALSE){
  labs=rownames(lts[[1]])
  nc=length(lts)
  nr=length(labs)
  if(set.seed.na){
    ssd<-NA
  } else {
    ssd<-runif(nr*nc,0,1)
  }
  pmat<-matrix(ssd,ncol=nc,nrow=nr,dimnames=list(labs,glabs))
  pmat<-as.data.frame(pmat)
  for(i in 1:length(lts)){
    l1<-lts[[i]][labs,"Pvalue"]
    idx<-!is.na(l1)
    pmat[idx,glabs[i]]<-l1[idx]
  }
  return(pmat)
}
##-----------------------------------------------------------------------------
#compute jaccard coefficient on tnets (among regulons)
getJC<-function(tnet){
  bin<-tnet
  bin[bin!=0]=1
  jc<-function(x,xmat){
    c<-x+xmat
    a<-colSums(c==2)
    b<-colSums(c>0)
    b[b==0]=1
    a/b
  }
  jcmat<-apply(bin,2,jc,xmat=bin)
  colnames(jcmat)<-colnames(bin)
  diag(jcmat)=0
  jcmat
}

##-----------------------------------------------------------------------------
##get adjacence matrices from dpi and a set of regulons
getmat<-function(reguls,dpimat,testsym=FALSE){
  ck1<-!reguls%in%colnames(dpimat)
  ck2<-!reguls%in%rownames(dpimat)
  if(sum(ck1)>0 || length(reguls)>ncol(dpimat)){
    stop("Olha direito isso!!!")
  }
  if(sum(ck2)>0){
    addr<-matrix(0,nrow=sum(ck2), ncol=ncol(dpimat))
    rownames(addr)<-reguls[ck2]
    dpimat<-rbind(addr,dpimat)
  }
  onlyreg<-reguls
  dpimat<-dpimat[,onlyreg]
  bin<-dpimat
  bin[bin!=0]=1
  idx<-apply(bin,1,sum)>0
  onlytar<-rownames(bin)[idx]
  onlytar<-onlytar[!onlytar%in%onlyreg]
  dpimat<-dpimat[c(onlyreg,onlytar),]
  tmat<-t(dpimat[onlytar,])
  zmat<-matrix(0,length(onlytar),length(onlytar))
  bmat<-rbind(tmat,zmat)
  adjmt<-cbind(dpimat,bmat)
  if(testsym){
    if(isSymmetric.matrix(adjmt,tol=0)){
      return(adjmt)
    } else {
      stop("result is not symmetric!!")
    }
  } else {
    return(adjmt)
  }
}
##-----------------------------------------------------------------------------
##color scale palette with breaks
colorscale1=function(x,breaks,cols,na.col,isrev,nquant,roundleg){
  # check arg
  if(!is.null(nquant)){
    if(nquant<2)stop("NOTE: require at least two quantiles!")
    breaks=quantile(x,probs=seq(0,1,length.out=nquant),na.rm=TRUE,names=FALSE)
  }
  if(is.null(cols))cols=c("darkblue","blue","orange","cyan","red","darkred")
  if(is.null(na.col)){na.col=grey(0.7)} else {na.col=na.col[1]}	
  if(is.character(x))stop("NOTE: 'breaks' arg. can not be applyed to characters!")
  if(sum(is.null(breaks))>0)stop("NOTE: breaks do not support null values!")
  if(length(breaks)<3)stop("NOTE: require at least three breaks!")
  if(isrev)cols=rev(cols)
  # adjust breaks and get palette
  bkcenter=(breaks[-length(breaks)]+breaks[-1])/2
  bkcenter=c(-Inf,bkcenter,+Inf)
  cols=colorRampPalette(colors=cols)(length(bkcenter)-1)
  # set colors to x
  x.col=rep(NA,length(x))
  cuts=cut(x[!is.na(x)], breaks=bkcenter,include.lowest=TRUE)
  x.col[!is.na(x)]=cols[as.integer(cuts)]
  x.col[is.na(x.col)]=colorRampPalette(colors=c(na.col,na.col))(1)
  # get intervals
  if(is.null(nquant)){
    interv=levels(cuts)
  } else {
    interv=seq(0,1,length.out=nquant+1)[-1]
    interv=paste(interv*100,"%",sep="")
  }		
  breaks=format(breaks, digits=roundleg, nsmall=roundleg)	
  leg=list(scale=cols,legend=breaks, interval=interv)
  res=list(res=x.col,leg=leg)	
  return(res)
}
##-----------------------------------------------------------------------------
#neg/pos color scale palette with breaks (left/right)
colorscale2=function(x,breaks,cols,na.col,isrev,nquant,roundleg, set2center=TRUE){
  # check args
  if(!is.null(nquant)){
    if(nquant<2)stop("NOTE: require at least two quantiles!")
    breaks=quantile(x,probs=seq(0,1,length.out=nquant),na.rm=TRUE,names=FALSE)
  }	
  if(is.null(cols))cols=c("darkblue","white","darkred")
  if(is.null(na.col)){na.col=grey(0.7)} else {na.col=na.col[1]}
  if(is.character(x))stop("NOTE: 'breaks' arg. can not be applyed to characters!")
  if(sum(is.null(breaks))>0)stop("NOTE: breaks do not support null values!")
  if(length(breaks)<3)stop("NOTE: require at least three breaks!")
  if(isrev)cols=rev(cols)
  if(set2center){
    bkcenter=(breaks[-length(breaks)]+breaks[-1])/2
  } else {
    bkcenter=breaks
  }	
  # check color vec		
  lt=length(cols)
  if(lt/2==as.integer(lt/2))lt=lt+1	
  cols=colorRampPalette(colors=cols)(lt)
  lfrt=as.integer(lt/2)+1
  # get neg/pos colors
  negCols=cols[1:lfrt]
  posCols=cols[lfrt:lt]
  ct.col=cols[lfrt]
  # check and adjust breaks
  lt=length(bkcenter)
  if(lt/2==as.integer(lt/2)){
    lf=lt/2
    rt=(lt/2)+1
    center=(bkcenter[lf]+bkcenter[rt])/2
    negBreaks=c(-Inf,bkcenter[1:lf],center)
    posBreaks=c(center,bkcenter[rt:lt],+Inf)
  } else {
    lfrt=as.integer(lt/2)+1
    center=bkcenter[lfrt]
    negBreaks=c(-Inf,bkcenter[1:lfrt])
    posBreaks=c(bkcenter[lfrt:lt],+Inf)
  }
  # set main palettes	
  negCols=colorRampPalette(colors=negCols)(length(negBreaks))[-length(negBreaks)]
  posCols=colorRampPalette(colors=posCols)(length(posBreaks))[-1]
  # set minor palettesscale
  na.col=colorRampPalette(colors=c(na.col,na.col))(1)
  ct.col=colorRampPalette(colors=c(ct.col,ct.col))(1)		
  # set colors to x
  x.col=rep(NA,length(x))
  idx=x<center & !is.na(x)
  negcuts=cut(x[idx],breaks=negBreaks,include.lowest=TRUE)
  x.col[idx]=negCols[as.integer(negcuts)]
  idx=x>center & !is.na(x)
  poscuts=cut(x[idx],breaks=posBreaks)
  x.col[idx]=posCols[as.integer(poscuts)]
  x.col[x==center]=ct.col
  x.col[is.na(x.col)]=na.col
  # get intervals
  if(is.null(nquant)){
    interv=c(levels(negcuts),levels(poscuts))
  } else {
    interv=seq(0,1,length.out=nquant+1)[-1]
    interv=paste(interv*100,"%",sep="")
  }
  testlen=length(breaks)/2
  if(as.integer(testlen)<testlen){
    idx=as.integer(testlen)+1
    breaks=breaks[c(1:idx,idx:length(breaks))]
  }
  breaks=c(format(breaks, digits=roundleg, nsmall=roundleg))
  leg=list(scale=c(negCols,posCols),legend=breaks, interval=interv)
  res=list(res=x.col,leg=leg)
  return(res)
}
##-----------------------------------------------------------------------------
##plot legends
plotleg1<-function(leg,title="", fname="leg"){
  xlab<-format(names(leg), digits = 3)
  pdf(file=paste(fname,".pdf",sep=""), width=5, height=2.5)
  image(matrix(1:length(leg)),col=leg,axes=FALSE)
  sq<-as.logical(1:length(xlab))
  axis(1, at=seq(0, 1, length.out = length(xlab))[sq], labels=xlab[sq], las=2,cex.axis=1)
  box()
  mtext(title,adj=0,line=1,cex=1.2)
  dev.off()
  cat(paste("Legend '", fname,".pdf' generated!\n\n", sep=""))
}
plotleg2<-function(leg, title="",cex=0.6){
  xlab<-format(leg$legend, digits=3)
  image(matrix(1:length(leg$scale)), col=leg$scale, axes=FALSE)
  sq<-!as.logical(1:length(xlab)%%2)
  axis(1, at=seq(0, 1, length.out=length(xlab))[sq], 
       labels=xlab[sq], las=2,cex.axis=cex,tcl=-0.2,lwd=0.5)
  box(lwd=0.5,col="grey50")
  mtext(title,adj=0,line=0.1,cex=0.5)
}
##-----------------------------------------------------------------------------
##plot overlap, synergy and shadow
plotovsysh<-function(ovmat,symat,shmat,xlab,leg1,leg2,fname="synergy.shadow"){ 
  fname=paste(fname,".pdf",sep="")
  pdf(file=fname,height=4.8,width=8.8) 
  layout(matrix(c(0,1,2,3,0,0,0,4,5,0), 2, 5, byrow = TRUE), heights=c(6,2), widths=c(1,3,3,3,1))
  par(mar=c(6, 1, 4, 0))
  #---plot overlap
  pmat<-matrix(ovmat,ncol=ncol(ovmat))
  image(matrix(1:prod(dim(pmat)),ncol=ncol(pmat)),col=pmat,axes=FALSE)
  axis(1, at=seq(0,1,length.out=ncol(pmat)), labels=xlab, las=2,
       cex.axis=1,lwd=0.5,tcl=-0.2,hadj=1,mgp=c(3,0.5,0))
  axis(2, at=seq(0,1,length.out=ncol(pmat)), labels=xlab, las=2,
       cex.axis=1,lwd=0.5,tcl=-0.2,hadj=1,mgp=c(3,0.5,0))
  box(lwd=0.75)
  mtext("Overlap",side=3, adj=0,line=0.5,cex=1.2)
  mtext("Regulons",side=2, adj=0.5,line=5,cex=0.9)
  mtext("Regulons",side=1, adj=0.5,line=5,cex=0.9)
  #---plot synergy
  pmat<-matrix(symat,ncol=ncol(symat))
  image(matrix(1:prod(dim(pmat)),ncol=ncol(pmat)),col=pmat,axes=FALSE)
  axis(1, at=seq(0,1,length.out=ncol(pmat)), labels=xlab, las=2,
       cex.axis=1,lwd=0.5,tcl=-0.2,hadj=1,mgp=c(3,0.5,0))
  box(lwd=0.75)
  mtext("Synergy",side=3, adj=0,line=0.5,cex=1.2)
  mtext("Regulons",side=1, adj=0.5,line=5,cex=0.9)
  #---plot shadow
  pmat<-matrix(shmat,ncol=ncol(shmat))
  image(matrix(1:prod(dim(pmat)),ncol=ncol(pmat)),col=pmat,axes=FALSE)
  abline(a=0, b=1, col = "black", lty=3,lwd=1.2)
  axis(1, at=seq(0,1,length.out=ncol(pmat)), labels=xlab, las=2,
       cex.axis=1,lwd=0.5,tcl=-0.2,hadj=1,mgp=c(3,0.5,0))
  box(lwd=0.75)
  mtext("Shadow",side=3, adj=0,line=0.5,cex=1.2)
  mtext("Shadows",side=1, adj=0.5,line=5,cex=0.9)
  #---plot legends
  par(mar = c(3.7, 2, 3, 2))
  #--Jaccard
  scl<-leg1$scale
  leg<-leg1$legend
  image(matrix(1:length(scl)),col=scl,axes=FALSE)
  sq<-as.logical(1:length(leg)%%2)
  axis(1, at=seq(0,1,length.out=length(leg))[sq], labels=leg[sq], las=2,
       cex.axis=1.0,lwd=0.5,tcl=-0.2,hadj=0.8)
  box(lwd=0.5)
  mtext("Jaccard coefficient",adj=0,line=0,cex=0.8)
  #--pvalue
  scl<-leg2$scale
  leg<-leg2$legend
  image(matrix(1:length(scl)),col=scl,axes=FALSE)
  sq<-as.logical(1:length(leg)%%2)
  axis(1, at=seq(0,1,length.out=length(leg))[sq], labels=leg[sq], las=2,
       cex.axis=1.0,lwd=0.5,tcl=-0.2,hadj=0.8)
  box(lwd=0.5)
  mtext("adj.pvalue",adj=0,line=0,cex=0.8)
  dev.off()
  cat(paste("File '", fname,"' generated!\n\n", sep=""))
}
##-----------------------------------------------------------------------------
##plot MRA summary as heatmap
plot.summary.mra<-function(tb){
  #---plot as table
  img<-t(tb[rev(1:nrow(tb)),])
  pdf(file="consensus.masters.pdf",height=4.2,width=3)
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(8,2))
  par(mar = c(0, 5, 6, 2), mgp=c(3,0.3,0))
  cols<-colorRampPalette(colors=brewer.pal(6,"Blues"))(100)[c(1,25,62,85)]
  image(img,col=cols,axes=FALSE)
  axis(2, at=seq(0,1,length.out=ncol(img)), labels=colnames(img), las=2,
       cex.axis=0.7,lwd=0.5,tcl=-0.2,hadj=1,mgp=c(3,0.5,0))
  axis(3, at=seq(0,1,length.out=nrow(img)), labels=rep(c("Exp1","Exp2","Exp3"),4), las=2,
       cex.axis=0.7,lwd=0.5,tcl=-0.2,hadj=0,mgp=c(3,0.5,0)) 
  box(lwd=0.5)
  mtext("Master regulators",side=2,adj=0.5,line=3.5,cex=0.7)
  mtext(text=c("","Cohort I","","","Cohort II","","","Normals","","","T-ALL",""), 
        at=seq(0,1,length.out=nrow(img)), side=3, cex=0.7, line=2.5)
  #---plot legends
  par(mar = c(2.5, 13.2, 1.7, 2),mgp=c(3,0.3,0))
  xlab<-c("0.00","0.01","0.05")
  image(matrix(1:length(xlab)), col=rev(cols)[1:3], axes=FALSE)
  sq<-as.logical(1:length(xlab))
  axis(1, at=seq(0, 1, length.out=length(xlab))[sq], 
       labels=xlab[sq], las=2,cex.axis=0.6,tcl=-0.2,lwd=0.5)
  box(lwd=0.5,col="grey50")
  mtext("DPI cutoff",adj=0,line=0.1,cex=0.6) 
  dev.off()
  cat("-File 'consensus.masters.pdf' generated!\n\n")
}
