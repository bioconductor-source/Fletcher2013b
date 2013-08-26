###############################################################################
# Motif analysis using bdsite maps derived from ChIP-seq data     			    
###############################################################################
getGeneDistBD<-function(tssmap, bdsites,chromlen){
  delta_bd <- matrix(NA,nrow=nrow(tssmap),ncol=5)
  colnames(delta_bd) <- c("nearup","neardown","nearbd","randbd","background")
  rownames(delta_bd) <- tssmap$entrezgene
  for(i in 1:nrow(tssmap)){
    chrom <- as.character(tssmap$chrom[[i]])
    if(chrom%in%names(chromlen)){
      txstart <- tssmap$transcript_start[i]
      strd <- ifelse(tssmap$strand[i]=="+",1,-1)
      #----
      bdstart<-bdsites$start[bdsites$chrom==chrom]
      bdend<-bdsites$end[bdsites$chrom==chrom]
      delta <- (bdstart - txstart) * strd
      delta1 <- (bdend - txstart) * strd
      idx<-abs(delta1)>abs(delta)
      delta[idx]<-delta1[idx]
      idx<-abs(delta+delta1)!=(abs(delta)+abs(delta1))
      delta[idx]<-0
      #----
      upstrdelta <- delta[delta<=0]
      dwstrdelta <- delta[delta>=0]
      delta_bd[i,1] <- ifelse(length(upstrdelta)>0,max(upstrdelta, na.rm=TRUE)/1000,NA)
      delta_bd[i,2] <- ifelse(length(dwstrdelta)>0, min(dwstrdelta, na.rm=TRUE)/1000,NA)
      idx<-sort.list(abs(delta))
      delta_bd[i,3] <- ifelse(length(delta)>0,delta[idx][1]/1000,NA)        
      #---
      bds <- runif(length(delta), min = 1, max = chromlen[[chrom]])
      #----
      len<-bdend-bdstart
      bdstart<-bds-len
      bdend<-bds+len
      delta <- (bdstart - txstart) * strd
      delta1 <- (bdend - txstart) * strd
      idx<-abs(delta1)>abs(delta)
      delta[idx]<-delta1[idx]
      idx<-abs(delta+delta1)!=(abs(delta)+abs(delta1))
      delta[idx]<-0
      idx<-sort.list(abs(delta))
      delta_bd[i,4] <- ifelse(length(delta)>0,delta[idx][1]/1000,NA)
      delta_bd[i,5] <- ifelse(length(delta)>0,delta[sample(idx)][1]/1000,NA)
    }
  }
  delta_bd <- data.frame(ENTREZ=rownames(delta_bd), 
                         nearup=delta_bd[,1], neardown=delta_bd[,2],
                         nearbd=delta_bd[,3],randbd=delta_bd[,4],
                         background=delta_bd[,5],stringsAsFactors=FALSE)
  rownames(delta_bd) <- delta_bd$ENTREZ
  delta_bd
}
##-----------------------------------------------------------------------------
getProbeDistBD<-function(geneids, delta_bd){
  probe.list <- data.frame(geneids, nearup=delta_bd[geneids$ENTREZ,"nearup"],
                           neardown=delta_bd[geneids$ENTREZ,"neardown"],nearbd=delta_bd[geneids$ENTREZ,"nearbd"],
                           randbd=delta_bd[geneids$ENTREZ,"randbd"],background=delta_bd[geneids$ENTREZ,"background"],
                           stringsAsFactors=FALSE)
  probe.list
}
##-----------------------------------------------------------------------------
getstat<-function(bin, bd.dist,nper=100){
  bdpvals<-rep(NA,ncol(bin))
  bdpvals.sites<-rep(NA,ncol(bin))
  names(bdpvals)<-colnames(bin)
  names(bdpvals.sites)<-colnames(bin)
  for(i in colnames(bin)){
    reg<-bin[,i]
    reg<-names(reg[reg==1])
    noreg<-rownames(bin)
    randmed<-rep(NA,nper)
    randmed.sites<-rep(NA,nper)
    for(j in 1:nper){
      rand<-noreg[sample(1:length(noreg))[1:length(reg)]]
      tp<-abs(bd.dist[rand,"nearbd"])
      randmed[j]<-median(tp,na.rm=TRUE)
      tp<-abs(bd.dist[rand,"randbd"])
      randmed.sites[j]<-median(tp,na.rm=TRUE)
    }
    tp<-abs(bd.dist[reg,"nearbd"])
    regmed<-median(tp,na.rm=TRUE)
    bdpvals[i]<-sum(regmed<randmed)
    bdpvals.sites[i]<-sum(regmed<randmed.sites)
  }
  bdpvals=(nper-bdpvals)/(nper)
  bdpvals.sites=(nper-bdpvals.sites)/(nper)
  list(pvreg=bdpvals,pvsite=bdpvals.sites)
}
##-----------------------------------------------------------------------------
getTSSmap<-function(ENTREZIDS){
  ids<-unique(ENTREZIDS)
  ids<-ids[!is.na(ids)]
  listMarts(host='may2009.archive.ensembl.org')
  ensembl=useMart(host='may2009.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL')
  listDatasets(ensembl)
  ensembl=useMart(host='may2009.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', 
                  dataset='hsapiens_gene_ensembl')
  mart <- useDataset("hsapiens_gene_ensembl", ensembl)
  tssmap <- getBM(attributes=c("entrezgene","chromosome_name","transcript_start","strand"), mart=mart)
  tssmap <- tssmap[tssmap$entrezgene%in%ids,]
  tssmap <- data.frame(tssmap, chrom=paste("chr", tssmap$chromosome_name,sep=""), stringsAsFactors=FALSE)
  tssmap[tssmap$strand=="1","strand"]="+"
  tssmap[tssmap$strand=="-1","strand"]="-"
  tp<-tssmap[,c(1,3)]
  tp<-aggregate(. ~ entrezgene, data=tp, mean, na.rm=TRUE)
  rownames(tp)<-tp$entrezgene
  idx<-match(rownames(tp),tssmap$entrezgene)
  tssmap<-data.frame(tp,tssmap[idx,c("chrom","strand")], stringsAsFactors=FALSE)
  tssmap
}
##-----------------------------------------------------------------------------
plotregs<-function(bd.dist, bdpvals, bin, allmasters, labmaster="ESR1", xlim=c(-400,400), 
                   ylim=c(0,0.01), nplot=200, file=NULL, adjust=1.2){
  idmaster<-allmasters[labmaster]
  pvreg <-format(bdpvals$pvreg, scientific = TRUE,digits = 2)
  pvsite <-format(bdpvals$pvsite, scientific = TRUE,digits = 2)
  pvreg<-pvreg[idmaster]
  pvsite<-pvsite[idmaster]
  plotrg<-function(bd.dist, idmaster, bin, labmaster, xlim, ylim, nplot, pvreg, pvsite, file, adjust){
    reg<-bin[,idmaster]
    reg<-names(reg[reg==1])	
    tp<-bd.dist[reg,]
    #---res0 - main results
    res0<-density(tp$nearbd,na.rm=TRUE, adjust=adjust, n=length(tp$nearbd),from=xlim[1],to=xlim[2])
    #---res1
    res1<-list()
    for(i in 1:nplot){
      randreg<-sample(rownames(bin))[1:length(reg)]
      tp<-bd.dist[randreg,]
      res1[[i]]<-density(tp$nearbd,na.rm=TRUE, adjust=adjust, n=length(tp$nearbd),from=xlim[1],to=xlim[2])
    }
    xr<-NULL
    yr<-NULL
    for(i in 1:nplot){
      xr<-cbind(xr,res1[[i]]$x)
      yr<-cbind(yr,res1[[i]]$y)
    }
    res1$x=apply(xr,1,mean)
    res1$y=apply(yr,1,mean)
    res1$sd=apply(yr,1,sd)
    #---res2
    res2<-list()
    for(i in 1:nplot){
      randreg<-sample(rownames(bin))[1:length(reg)]
      tp<-bd.dist[randreg,]
      res2[[i]]<-density(tp$randbd,na.rm=TRUE, adjust=adjust, n=length(tp$randbd),from=xlim[1],to=xlim[2])
    }
    xr<-NULL
    yr<-NULL
    for(i in 1:nplot){
      xr<-cbind(xr,res2[[i]]$x)
      yr<-cbind(yr,res2[[i]]$y)
    }
    res2$x=apply(xr,1,mean)
    res2$y=apply(yr,1,mean)
    res2$sd=apply(yr,1,sd)   
    #---res3 - background
    res3<-density(tp$background,na.rm=TRUE, adjust=adjust, n=length(tp$background),from=xlim[1],to=xlim[2])
    #---PLOTS---
    if(!is.null(file)){
      pdf(file=paste(file,".pdf",sep=""), height=3.0, width=2.4)
    }
    plot.new()
    plot.window(ylim=ylim,xlim=xlim)
    axis(2,cex.axis=0.5,las=1,tcl=-0.3) 
    axis(1,cex.axis=0.5,las=1,tcl=-0.3)
    title(xlab="Nearest binding site (kb)",ylab="TSS density",cex.lab=0.8) 
    #---plot3
    lines(x=res3$x, y=res3$y, xlim=xlim,col="gray60", lwd=1)
    #---areas
    polygon(c(res1$x,rev(res1$x)), c(res1$y+res1$sd,rev(res1$y-res1$sd)), col="gray90", border="gray85",lwd=0.3)
    polygon(c(res2$x,rev(res2$x)), c(res2$y+res2$sd,rev(res2$y-res2$sd)), col="gray90", border="gray85",lwd=0.3)    
    #---plot1
    #lines(x=res1$x, y=res1$y+res1$sd, xlim=xlim,col="palegreen3", lwd=0.2)	
    lines(x=res1$x, y=res1$y, xlim=xlim,col="green2", lwd=0.8)
    #lines(x=res1$x, y=res1$y-res1$sd, xlim=xlim,col="palegreen3", lwd=0.2) 
    #---plot2
    #lines(x=res2$x, y=res2$y+res2$sd, xlim=xlim,col="black", lwd=0.2)	
    lines(x=res2$x, y=res2$y, xlim=xlim,col="blue1", lwd=0.8)
    #lines(x=res2$x, y=res2$y-res2$sd, xlim=xlim,col="black", lwd=0.2)
    #---plot0 (main results)
    lines(x=res0$x, y=res0$y, xlim=xlim,col="red", lwd=1)
    #---legends
    lab1<-paste(labmaster," regulon",sep="")
    lab2<-paste("random regulons (p=", pvreg,")",sep="")
    if(as.numeric(pvreg)==0)
      lab2<-paste("random regulons (p<", format(1/1000, scientific = TRUE,digits = 2),")",sep="")
    lab3<-paste("random sites (p=", pvsite,")",sep="")
    if(as.numeric(pvsite)==0)
      lab3<-paste("random sites (p<", format(1/1000, scientific = TRUE,digits = 2),")",sep="")
    legend("topleft",c(lab1,lab2,lab3,"background"),bty="n",cex=0.35,
           col=c("red", "green2","blue1","gray60"), pch="---", pt.cex=0.8)
    if(!is.null(file)){
      dev.off()
      cat(paste("File '", paste(file,".pdf",sep=""),"' generated!\n\n", sep=""))
    }
  }
  plotrg(bd.dist, idmaster, bin, labmaster, xlim, ylim, nplot, pvreg, pvsite, file, adjust)
}
##-----------------------------------------------------------------------------
plotallbd<-function(dist.list,reg,labmaster, pvreg, xlim=c(-400,400), ylim=NULL, file=NULL){
  xr<-NULL
  yr<-NULL    
  pvreg <-format(pvreg, scientific = TRUE,digits = 2)
  for(i in 1:200){
    tp<-sample(rownames(dist.list[["rdmotifs"]]))[1:length(reg)]
    tp<-dist.list[["rdmotifs"]][tp,]
    tp<-density(tp$nearbd,na.rm=TRUE, adjust=1.5, n=nrow(tp),from=xlim[1],to=xlim[2]) 
    xr<-cbind(xr,tp$x)
    yr<-cbind(yr,tp$y)      
  }
  xmd=apply(xr,1,mean)
  ymd=apply(yr,1,mean)
  ysd=apply(yr,1,sd)  
  #------------------
  res<-list()
  for(i in names(dist.list)){
    tp<-dist.list[[i]][reg,]    
    if(sum(!is.na(tp$nearbd))>3){
      nr=length(tp$nearbd)
      res[[i]]<-density(tp$nearbd,na.rm=TRUE, adjust=1.5, n=nr,from=xlim[1],to=xlim[2])      
    } else {
      res[[i]]$x<-xlim
      res[[i]]$y<-c(0,0) 
    }
  }
  #------------------
  maxy<-max(c(res[[1]]$y, res[[2]]$y),na.rm=TRUE)*1.5
  if(is.null(ylim))ylim=c(0,maxy)
  if(!is.null(file))pdf(file=paste(file,".pdf",sep=""), height=3.0, width=2.4)
  plot.new()
  plot.window(ylim=ylim,xlim=xlim)
  axis(2,cex.axis=0.5,las=1,tcl=-0.3) 
  axis(1,cex.axis=0.5,las=1,tcl=-0.3)
  title(xlab="Nearest binding site (kb)",ylab="TSS density",cex.lab=0.8)
  #lines(x=xmd, y=ymd+ysd, xlim=xlim,col="grey80", lwd=0.7)
  polygon(c(xmd,rev(xmd)), c(ymd+ysd,rev(ymd-ysd)), col="gray90", border="gray85",lwd=0.3)
  #lines(x=res[[2]]$x, y=res[[2]]$y, xlim=xlim,col="black", lwd=1.5)
  lines(x=xmd, y=ymd, xlim=xlim,col="black", lwd=1.2)
  #lines(x=xmd, y=ymd-ysd, xlim=xlim,col="grey80", lwd=0.7)
  lines(x=res[[1]]$x, y=res[[1]]$y, xlim=xlim,col="red", lwd=1.2)
  #legends
  lab1<-paste(labmaster," regulon",sep="")
  lab2<-paste("random PWM (p=", pvreg,")",sep="")
  legend("topleft",c(lab1,lab2),bty="n",cex=0.35,col=c("red", "black"), pch="---", pt.cex=0.8)
  if(!is.null(file)){
    dev.off()
    cat(paste("File '", paste(file,".pdf",sep=""),"' generated!\n\n", sep=""))
  }
}
