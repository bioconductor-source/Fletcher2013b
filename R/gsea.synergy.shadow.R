###############################################################################
#Compute synergy and shadow analysis for Exp1              
###############################################################################
Fletcher2013pipeline.synergyShadow<-function(){
  cat("Running analysis pipeline ... \n\n")
  data(rtni1st,rtniIDs)
  data(miscellaneous)
  rtni1st<-get("rtni1st")
  rtniIDs<-get("rtniIDs")
  consensus<-get("consensus")
  masters<-consensus$masters
  data(Exp1limma)
  dt<-get("Exp1limma")
  #get pheno e hits for Exp1
  pheno<-dt$deg[,"coef.t6.E2FGF10-t6.E2"]
  names(pheno)<-rownames(dt$deg)
  phenoIDs<-dt$deg[,c("PROBEID","ENTREZ")]
  #run synergy analysis
  names(pheno)<-phenoIDs[names(pheno),"ENTREZ"]
  pheno<-pheno[!is.na(names(pheno))]
  rtni1st@results$tn.dpi<-rtni1st@results$tn.dpi[,masters]
  rtni1st@results$tn.ref<-rtni1st@results$tn.ref[,masters]
  rtni1st@transcriptionFactors=masters
  rtna <- tni2tna.preprocess(rtni1st, phenotype=pheno, verbose=FALSE)
  
  rtna <- tna.overlap(rtna, pValueCutoff=0.001, minRegulonSize=20, verbose=TRUE)
  rtna <- tna.synergy(rtna, pValueCutoff=0.001, minRegulonSize=20, nPermutations=1000, stepFilter=FALSE, verbose=TRUE)
  rtna <- tna.shadow(rtna, pValueCutoff=0.001, minRegulonSize=20, nPermutations=1000, stepFilter=FALSE, verbose=TRUE)
  plotSynergyShadow(rtna)
}
##-----------------------------------------------------------------------------
plotSynergyShadow<-function(rtna){
  data(rtni1st)
  data(miscellaneous)
  rtni1st<-get("rtni1st")
  consensus<-get("consensus")
  masters<-consensus$masters
  #compute jaccard coefficient and clustering
  jcmast<-abs(rtni1st@results$tn.ref[,masters])
  masters[]<-rtni1st@annotation[masters,"ENTREZ"]
  colnames(jcmast)<-masters
  jcmast<-getJC(jcmast)
  diag(jcmast)<-1
  hc<-hclust(as.dist(1-cor(jcmast)), method='ward')
  dend <- as.dendrogram(hc)
  dend<-reorder(x=dend, wts=rowMeans(jcmast),agglo.FUN=max)
  jcmast<-jcmast[order.dendrogram(dend),order.dendrogram(dend)]
  masters<-masters[order.dendrogram(dend)]
  #set commom colors and breaks
  #pvalues
  cpall<-c(colorRampPalette(c("darkred", "orange","grey98"),bias=2)(16)[c(2:6,15)],colorRampPalette(c("grey98", "grey80"))(7))
  bksPvals<-sort(c(10^(-rev(c(1:7))),c(2:5)*0.1,1,0.05))
  #overlap
  coall<-colorRampPalette(c("grey98",brewer.pal(6,"Blues")),bias=1.5)(11)
  bksOverlap<-seq(0.0,0.5,length.out=11)
  #get overlap matrix
  ovmat<-jcmast
  ovleg<-colorscale1(ovmat, breaks=bksOverlap, nquant=NULL, cols=coall, na.col="transparent", isrev=FALSE, roundleg=2)
  ovmat[,]<-ovleg$res
  #get synergy matrix  
  syn<-rtna@results$synergy.results
  symat<-jcmast
  symat[,]=1
  sapply(rownames(symat),function(i){
    idx<-which(syn[,"Regulon1"]==i)
    j<-syn[idx,"Regulon2"]
    symat[i,j]<<-syn[idx,"Adjusted.Pvalue"]
    symat[j,i]<<-syn[idx,"Adjusted.Pvalue"]
    NULL
  })
  syleg<-colorscale2(symat, breaks=bksPvals, nquant=NULL, cols=cpall, na.col="transparent", isrev=FALSE, roundleg=2, set2center=FALSE)
  symat[,]<-syleg$res
  #get shadow matrix  
  #---R1
  shad<-rtna@results$shadow.results$shadowR1
  shmat1<-jcmast
  shmat1[,]=1
  sapply(rownames(shmat1),function(i){
    idx<-which(shad[,"Regulon"]==i)
    j<-shad[idx,"Shadow"]
    shmat1[i,j]<<-shad[idx,"Adjusted.Pvalue"]
    shmat1[j,i]<<-shad[idx,"Adjusted.Pvalue"]
  })
  shmat1[lower.tri(shmat1)]<-1
  #---R2
  shad<-rtna@results$shadow.results$shadowR2
  shmat2<-jcmast
  shmat2[,]=1
  sapply(rownames(shmat2),function(i){
    idx<-which(shad[,"Regulon"]==i)
    j<-shad[idx,"Shadow"]
    shmat2[i,j]<<-shad[idx,"Adjusted.Pvalue"]
    shmat2[j,i]<<-shad[idx,"Adjusted.Pvalue"]
  })
  shmat2[upper.tri(shmat2)]<-1
  shmat<-shmat1*shmat2
  shleg<-colorscale2(shmat, breaks=bksPvals, nquant=NULL, cols=cpall, na.col="transparent", isrev=FALSE, roundleg=2, set2center=FALSE)
  shmat[,]<-shleg$res
  #Plot overlap,synergy and shadow
  plotovsysh(ovmat,symat,shmat,names(masters),ovleg$leg,syleg$leg)
}
##-----------------------------------------------------------------------------
Fletcher2013gsea.regulons<-function(what="Exp1", timepoint=6, verbose=TRUE){
  cat("Running analysis pipeline ... \n\n")
  if(sum(what%in%c("Exp1","Exp2","Exp3"))!=1)stop("'what' should be any of 'Exp1','Exp2' and 'Exp3'!")
  if(sum(timepoint%in%c(6,24))!=1)stop("'timepoint' should 6 or 24!")
  data(rtni1st,rtniIDs)
  data(miscellaneous)
  rtni1st<-get("rtni1st")
  rtniIDs<-get("rtniIDs")
  consensus<-get("consensus")
  masters<-consensus$masters
  consensusmasters<-consensus$consensusmasters
  #get pheno e hits for Exp1
  if(what=="Exp1"){
    data(Exp1limma)
    dt<-get("Exp1limma")
    if(timepoint==24){
      pheno<-dt$deg[,"coef.t24.E2FGF10-t24.E2"]
      labPheno="E2.FGF10.24h"
    } else {
      pheno<-dt$deg[,"coef.t6.E2FGF10-t6.E2"]
      labPheno="E2.FGF10.6h"
    }
    names(pheno)<-rownames(dt$deg)
    phenoIDs<-dt$deg[,c("PROBEID","ENTREZ")]
    tnetIDs=rtniIDs[,c("PROBEID","ENTREZ")]
  } else if(what=="Exp2"){
    data(Exp2limma)
    dt<-get("Exp2limma")
    if(timepoint==24){
      pheno<-dt$deg[,"coef.t24.E2.AP20187-t24.E2"]
      labPheno="E2.AP20187.24h"
    } else {
      pheno<-dt$deg[,"coef.t6.E2.AP20187-t6.E2"]
      labPheno="E2.AP20187.6h"
    }
    names(pheno)<-rownames(dt$deg)
    phenoIDs<-dt$deg[,c("PROBEID","ENTREZ")]
    tnetIDs=rtniIDs[,c("PROBEID","ENTREZ")] 
  } else if(what=="Exp3"){
    data(Exp3limma)
    dt<-get("Exp3limma")
    if(timepoint==24){
      pheno<-dt$deg[,"coef.PlusTet.PlusE2PlusFGF10.24h-PlusTet.PlusE2.24h"]
      labPheno="Tet.E2.FGF10.24h"
    } else {
      pheno<-dt$deg[,"coef.PlusTet.PlusE2PlusFGF10.6h-PlusTet.PlusE2.6h"]
      labPheno="Tet.E2.FGF10.6h"
    }
    names(pheno)<-rownames(dt$deg)
    phenoIDs<-dt$deg[,c("PROBEID","ENTREZ")]
    tnetIDs=rtniIDs[,c("PROBEID","ENTREZ")]
  }
  #run gsea1 analysis
  names(pheno)<-phenoIDs[names(pheno),"ENTREZ"]
  pheno<-pheno[!is.na(names(pheno))]
  rtni1st@results$tn.dpi<-rtni1st@results$tn.dpi[,consensusmasters]
  rtni1st@results$tn.ref<-rtni1st@results$tn.ref[,consensusmasters]
  rtni1st@transcriptionFactors=consensusmasters
  rtna <- tni2tna.preprocess(rtni1st, phenotype=pheno, verbose=FALSE)
  
  rtna <- tna.gsea1(rtna, nPermutations=1000, stepFilter=FALSE, verbose=verbose)
  #plot gsea1
  palette(c("black","cyan","red","green3","blue"))
  tna.plot.gsea1(rtna,ylimPanels=c(0,3.5,0,0.6),labPheno=labPheno, ntop=-1)
  cat(paste("File '", labPheno,".pdf' generated!\n\n", sep=""))
}

