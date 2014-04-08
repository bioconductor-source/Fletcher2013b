###############################################################################
#PLot the enrichmen map in an heatmap        	    
###############################################################################
Fletcher2013mra.heatmap1<-function(){
  cat("Running analysis pipeline ... \n\n")
  data(rtni1st,rtniIDs)
  rtni1st<-get("rtni1st")
  rtniIDs<-get("rtniIDs")
  #---compute jaccard coefficient
  minRegulonSize=20
  tn.ref<-abs(rtni1st@results$tn.ref)
  tn.ref<-tn.ref[,colSums(tn.ref!=0)>=minRegulonSize]
  jcmat<-getJC(tn.ref)
  #--- get signatures
  Exp1<-Fletcher2013pipeline.deg(what="Exp1")
  Exp2<-Fletcher2013pipeline.deg(what="Exp2")
  Exp3<-Fletcher2013pipeline.deg(what="Exp3")
  #--- get MRA
  resExp1<-list()
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp1$E2FGF10, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp1$E2FGF10 <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp1$E2, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp1$E2 <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)  
  resExp2<-list()
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp2$E2AP20187, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp2$E2FGF10 <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp2$E2, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp2$E2 <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)  
  resExp3<-list()
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp3$TetE2FGF10, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp3$E2FGF10 <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp3$TetE2, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp3$E2 <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)  
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp3$random, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp3$random <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)  
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp3$random, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp3$random <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)
  #--- get signature aligned with jcmat
  getsigt<-function(jcmat,resExp1,resExp2,resExp3,bksPvals,cpall){
    #--- get name vector aline with jcmat
    alinedvec<-rep(NA,ncol(jcmat))
    names(alinedvec)<-colnames(jcmat)   
    #--- FGFR2 signature pvals  
    colsFGF10 <-list()
    colsFGF10[["Exp1 "]]<-alinedvec
    colsFGF10[["Exp1 "]][rownames(resExp1$E2FGF10)]<-resExp1$E2FGF10[,"Adjusted.Pvalue"]
    colsFGF10[["Exp2 "]]<-alinedvec
    colsFGF10[["Exp2 "]][rownames(resExp2$E2FGF10)]<-resExp2$E2FGF10[,"Adjusted.Pvalue"]  
    colsFGF10[["Exp3 "]]<-alinedvec
    colsFGF10[["Exp3 "]][rownames(resExp3$E2FGF10)]<-resExp3$E2FGF10[,"Adjusted.Pvalue"]
    colsFGF10[["random "]]<-alinedvec
    colsFGF10[["random "]][rownames(resExp3$random)]<-resExp3$random[,"Adjusted.Pvalue"]
    #--- E2 signature pvals
    colsE2<-list()
    colsE2[["E2 (Ctr1) "]]<-alinedvec
    colsE2[["E2 (Ctr1) "]][rownames(resExp1$E2)]<-resExp1$E2[,"Adjusted.Pvalue"]
    colsE2[["E2 (Ctr2) "]]<-alinedvec
    colsE2[["E2 (Ctr2) "]][rownames(resExp2$E2)]<-resExp2$E2[,"Adjusted.Pvalue"]  
    colsE2[["E2 (Ctr3) "]]<-alinedvec
    colsE2[["E2 (Ctr3) "]][rownames(resExp3$E2)]<-resExp3$E2[,"Adjusted.Pvalue"]
    colsE2[["random "]]<-alinedvec
    colsE2[["random "]][rownames(resExp3$random)]<-resExp3$random[,"Adjusted.Pvalue"]
    #--- set colors and breaks
    for(i in names(colsFGF10)){
      colsFGF10[[i]]<-colorscale2(colsFGF10[[i]], breaks=bksPvals, nquant=NULL, cols=cpall, 
                                  na.col=cpall[length(cpall)], isrev=FALSE, roundleg=2)$res
    }
    for(i in names(colsE2)){
      colsE2[[i]]<-colorscale2(colsE2[[i]], breaks=bksPvals, nquant=NULL, cols=cpall, 
                               na.col=cpall[length(cpall)], isrev=FALSE, roundleg=2)$res
    }
    list(colsFGF10=colsFGF10,colsE2=colsE2)
  }
  #set commom colors and breaks
  #pvalues
  cpall<-c(colorRampPalette(c("darkred", "orange","grey98"),bias=2)(16)[c(2:6,10)],colorRampPalette(c("grey98", "grey80"))(7))
  bksPvals<-sort(c(0,10^(-rev(c(1:6*2))),c(1:5)*0.1,1))
  #overlap
  coall<-colorRampPalette(c("grey98",brewer.pal(6,"Blues")),bias=1.5)(11)
  bksOverlap<-seq(0.0,0.5,length.out=11)
  #plot all regulons
  #get signatures aligned with jcmat (as side colors)
  sigt<-getsigt(jcmat,resExp1,resExp2,resExp3,bksPvals,cpall)  
  ids<-rtniIDs[colnames(jcmat),3]
  #set legends
  legends=list()
  leg<-colorscale2(bksPvals, breaks=bksPvals, nquant=NULL, cols=cpall, na.col="white", isrev=FALSE, roundleg=2)$leg
  legends[["Enrichment pvals"]]=leg$scale
  names(legends[["Enrichment pvals"]])<-leg$legend
  #set legends
  legends[['Jaccard coef.']]=coall
  names(legends[["Jaccard coef."]])<-bksOverlap
  #plot 1
  resol<-3.5 #set scale of plot size/resolution
  png(filename="heatmap_enrichmentmaps.png", width=500*resol, height=480*resol, res = 72*resol)
  dot4heatmap(jcmat, symm=TRUE,colorVec=coall, plotType="points", LegSideColors=legends,
              ColSideColors=sigt$colsFGF10, RowSideColors=sigt$colsE2,labCol=ids, 
              breaks=bksOverlap,labRowCol=ids, distfun=function(x) dist(x,method="manhattan"), plotScale=TRUE, 
              xlab="Regulons (network overlap)", ylab="", cexRowCol=1.2, pointSz=1.5,
              reorderfun=function(d,w) reorder(d, w), hclustfun=function(x) hclust(x, method='ward'),
              bgplot=coall[length(coall)],cexLeg=0.75, cexLab=1.4,sidewidth=0.11,margins=c(8,8),
              legwidth=0.32,matwidth=3.0)
  dev.off()
  cat("-File 'heatmap_enrichmentmaps.png' generated!\n\n")
  #simplify
  simplify<-function(){
    sigt<-getsigt(jcmat,resExp1,resExp2,resExp3,bksPvals,cpall)  
    ids<-rtniIDs[colnames(jcmat),3]
    resol<-3 #set scale of plot size/resolution
    names(sigt$colsE2)=""
    names(sigt$colsFGF10)=""
    png(filename="heatmap_enrichmentmaps.png", width=500*resol, height=480*resol, res = 72*resol)
    dot4heatmap(jcmat, symm=TRUE,colorVec=coall, plotType="points", ColSideColors=sigt$colsFGF10, 
                RowSideColors=sigt$colsE2, labCol="", cexLeg=0.85, cexRowCol=1.6,
                breaks=bksOverlap,labRowCol="", distfun=function(x) dist(x,method="manhattan"), plotScale=FALSE, 
                xlab="", ylab="", hclustfun=function(x) hclust(x, method='ward'),
                reorderfun=function(d,w) reorder(d, w,max), pointSz=1.5,
                sidewidth=0.11, dendwidth=0.3, margins=c(1,1),legwidth=0.32, matwidth=3.0)
    dev.off()
    cat("-File 'heatmap_enrichmentmaps.png' generated!\n\n")
  }
  #simplify()
  #legends
  #plotleg1(legends[[1]],"adj.pvalue","leg.enrichment")
  #plotleg1(legends[[2]],"Jaccard index","leg.overlap")
}
##-----------------------------------------------------------------------------
Fletcher2013mra.heatmap2<-function(){
  cat("Running analysis pipeline ... \n\n")
  data(rtni1st)
  data(miscellaneous)
  rtni1st<-get("rtni1st")
  consensus<-get("consensus")
  consensusmasters<-consensus$consensusmasters
  minRegulonSize=20
  #compute jaccard coefficient
  tn.ref<-abs(rtni1st@results$tn.ref[,consensusmasters])
  jcmast<-getJC(tn.ref)
  #get signatures
  Exp1<-Fletcher2013pipeline.deg(what="Exp1")
  Exp2<-Fletcher2013pipeline.deg(what="Exp2")
  Exp3<-Fletcher2013pipeline.deg(what="Exp3")
  #get MRA
  resExp1<-list()
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp1$E2FGF10, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp1$E2FGF10 <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp1$E2, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp1$E2 <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)  
  resExp2<-list()
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp2$E2AP20187, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp2$E2FGF10 <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp2$E2, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp2$E2 <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)  
  resExp3<-list()
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp3$TetE2FGF10, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp3$E2FGF10 <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp3$TetE2, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp3$E2 <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)  
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp3$random, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp3$random <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)  
  rtna <- tni2tna.preprocess(rtni1st, hits=Exp3$random, verbose=FALSE)
  rtna <- tna.mra(rtna, minRegulonSize=minRegulonSize, pAdjustMethod="holm", tnet="ref", verbose=FALSE)
  resExp3$random <- tna.get(rtna,what="mra", ntop=-1,reportNames=FALSE)
  #get signature aligned with jc
  getsigt<-function(jcmat,resExp1,resExp2,resExp3,bksPvals,cpall){
    #get name vector aline with jc mat
    alinedvec<-rep(NA,ncol(jcmat))
    names(alinedvec)<-colnames(jcmat)   
    #FGFR2 signature pvals  
    colsFGF10 <-list()
    colsFGF10[["Exp1 "]]<-alinedvec
    colsFGF10[["Exp1 "]][rownames(resExp1$E2FGF10)]<-resExp1$E2FGF10[,"Adjusted.Pvalue"]
    colsFGF10[["Exp2 "]]<-alinedvec
    colsFGF10[["Exp2 "]][rownames(resExp2$E2FGF10)]<-resExp2$E2FGF10[,"Adjusted.Pvalue"]  
    colsFGF10[["Exp3 "]]<-alinedvec
    colsFGF10[["Exp3 "]][rownames(resExp3$E2FGF10)]<-resExp3$E2FGF10[,"Adjusted.Pvalue"]
    colsFGF10[["random "]]<-alinedvec
    colsFGF10[["random "]][rownames(resExp3$random)]<-resExp3$random[,"Adjusted.Pvalue"]
    #E2 signature pvals
    colsE2<-list()
    colsE2[["E2 (Ctr1) "]]<-alinedvec
    colsE2[["E2 (Ctr1) "]][rownames(resExp1$E2)]<-resExp1$E2[,"Adjusted.Pvalue"]
    colsE2[["E2 (Ctr2) "]]<-alinedvec
    colsE2[["E2 (Ctr2) "]][rownames(resExp2$E2)]<-resExp2$E2[,"Adjusted.Pvalue"]  
    colsE2[["E2 (Ctr3) "]]<-alinedvec
    colsE2[["E2 (Ctr3) "]][rownames(resExp3$E2)]<-resExp3$E2[,"Adjusted.Pvalue"]
    colsE2[["random "]]<-alinedvec
    colsE2[["random "]][rownames(resExp3$random)]<-resExp3$random[,"Adjusted.Pvalue"]
    #set colors and breaks
    for(i in names(colsFGF10)){
      colsFGF10[[i]]<-colorscale2(colsFGF10[[i]], breaks=bksPvals, nquant=NULL, cols=cpall, 
                                  na.col=cpall[length(cpall)], isrev=FALSE, roundleg=2)$res
    }
    for(i in names(colsE2)){
      colsE2[[i]]<-colorscale2(colsE2[[i]], breaks=bksPvals, nquant=NULL, cols=cpall, 
                               na.col=cpall[length(cpall)], isrev=FALSE, roundleg=2)$res
    }
    list(colsFGF10=colsFGF10,colsE2=colsE2)
  }
  #set commom colors and breaks
  #pvalues
  cpall<-c(colorRampPalette(c("darkred", "orange","grey98"),bias=2)(16)[c(2:6,10)],colorRampPalette(c("grey98", "grey80"))(7))
  bksPvals<-sort(c(0,10^(-rev(c(1:6*2))),c(1:5)*0.1,1))
  #overlap
  coall<-colorRampPalette(c("grey98",brewer.pal(6,"Blues")),bias=1.5)(11)
  bksOverlap<-seq(0.0,1.0,length.out=11)
  #plot only consensusmasters
  diag(jcmast)<-1
  sigt<-getsigt(jcmast,resExp1,resExp2,resExp3,bksPvals,cpall)
  #set legends
  legends=list()
  leg<-colorscale2(bksPvals, breaks=bksPvals, nquant=NULL, cols=cpall, na.col="white", isrev=FALSE, roundleg=2)$leg
  legends[["Enrichment pvals"]]=leg$scale
  names(legends[["Enrichment pvals"]])<-leg$legend 
  #set legends
  legends[['Jaccard coef.']]=coall
  names(legends[["Jaccard coef."]])<-bksOverlap
  #test overlap
  #ph<-Fletcher2013.overlap(consensusmasters)[,1:2]
  #simplify
  simplify<-function(){
    hc<-hclust(dist(jcmast), method='ward')
    pdf(file="heatmap_enrichmentmap.master.pdf", width=4, height=3.8)
    dot4heatmap(x=jcmast,px=NULL, symm=TRUE,colorVec=coall, plotType="standard",
                ColSideColors=sigt$colsFGF10, RowSideColors=sigt$colsE2,labCol=names(consensusmasters), 
                cexRowCol=1.4, breaks=bksOverlap,labRowCol=names(consensusmasters),sidewidth=0.15,
                distfun=function(x) as.dist(1-cor(x)), plotScale=FALSE, ylab="",dendwidth=0.3,
                xlab="Enriched regulons (overlap)",reorderfun=function(d,w) reorder(d, w, max),
                hclustfun=function(x) hclust(x, method='ward'), matwidth=3.0,legwidth=0.32, 
                margins=c(7,7), cexLeg=0.75, cexLab=1.2, comment="")
    dev.off()
    cat("-File 'heatmap_enrichmentmap.master.pdf' generated!\n\n")
  }
  simplify()
  #legends
  #plotleg1(legends[[1]],"adj.pvalue","leg.enrichment")
  #plotleg1(legends[[2]],"Jaccard index","leg.overlap")
}

