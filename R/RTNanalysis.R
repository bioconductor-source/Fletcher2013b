###############################################################################
#MRA pipelines
###############################################################################
Fletcher2013pipeline.mra1st<-function(
  hits, minRegulonSize=20, idtype="probeid", pAdjustMethod="holm",
  tnet="dpi", eps=0, pValueCutoff=1e-4, verbose=TRUE, ...){
  cat("Running analysis pipeline ... \n\n")
  data(rtni1st)
  data(rtniIDs)
  rtni1st<-get("rtni1st")
  rtniIDs<-get("rtniIDs") 
  if(sum(tnet%in%c("dpi","ref"))!=1)stop("'tnet' should be any of 'dpi' and 'ref'!")
  if(sum(idtype%in%c("entrez","probeid"))!=1)stop("'idtype' should be any of 'probeid' and 'entrez'!")
  if(eps>0)rtni1st<-tni.dpi.filter(rtni1st, eps=eps)
  if(idtype=="entrez"){
    tnetIDs=rtniIDs[,1:2]
  } else {
    tnetIDs=NULL
  }
  hits<-hits[!is.na(hits)]
  rtna <- tni2tna.preprocess(rtni1st, hits=hits, verbose=FALSE)
  #compute regulons from both dpi and reference networks
  regs <- tna.mra(rtna, pValueCutoff=pValueCutoff, minRegulonSize=minRegulonSize, 
                   pAdjustMethod=pAdjustMethod, tnet=tnet, verbose=verbose)
  regs <- tna.get(regs,what="mra", ...=...)
  #remove regulons not listed in regref
  if(tnet=="dpi"){
    regref <- tna.mra(rtna, pValueCutoff=pValueCutoff, minRegulonSize=minRegulonSize, 
                     pAdjustMethod=pAdjustMethod, tnet="ref", verbose=verbose)
    regref <- tna.get(regref,what="mra", ...=...)
    idx<-rownames(regs)%in%rownames(regref)
    return(regs[idx,])
  } else {
    return(regs)
  }
}
##-----------------------------------------------------------------------------
Fletcher2013pipeline.mra2nd<-function(
  hits, minRegulonSize=20, idtype="probeid", pAdjustMethod="holm",
  tnet="dpi", eps=0, pValueCutoff=1e-4, verbose=TRUE, ...){
  cat("Running analysis pipeline ... \n\n")
  data(rtni2nd)
  data(rtniIDs)
  rtni2nd<-get("rtni2nd")
  rtniIDs<-get("rtniIDs")
  if(sum(tnet%in%c("dpi","ref"))!=1)stop("'tnet' should be any of 'dpi' and 'ref'!")
  if(sum(idtype%in%c("entrez","probeid"))!=1)stop("'idtype' should be any of 'probeid' and 'entrez'!")
  if(eps>0)rtni2nd<-tni.dpi.filter(rtni2nd, eps=eps)
  if(idtype=="entrez"){
    tnetIDs=rtniIDs[,1:2]
  } else {
    tnetIDs=NULL
  }
  hits<-hits[!is.na(hits)]
  rtna <- tni2tna.preprocess(rtni2nd, hits=hits, verbose=FALSE)
  #compute regulons from both dpi and reference networks
  regs <- tna.mra(rtna, pValueCutoff=pValueCutoff, minRegulonSize=minRegulonSize, 
                  pAdjustMethod=pAdjustMethod, tnet=tnet, verbose=verbose)
  regs <- tna.get(regs,what="mra", ...=...)
  #remove regulons not listed in regref
  if(tnet=="dpi"){
    regref <- tna.mra(rtna, pValueCutoff=pValueCutoff, minRegulonSize=minRegulonSize, 
                      pAdjustMethod=pAdjustMethod, tnet="ref", verbose=verbose)
    regref <- tna.get(regref,what="mra", ...=...)
    idx<-rownames(regs)%in%rownames(regref)
    return(regs[idx,])
  } else {
    return(regs)
  }
}
##-----------------------------------------------------------------------------
Fletcher2013pipeline.mraNormals<-function(
  hits, minRegulonSize=20, idtype="probeid", pAdjustMethod="holm",
  tnet="dpi", eps=0, pValueCutoff=1e-4, verbose=TRUE, ...){
  cat("Running analysis pipeline ... \n\n")
  data(rtniNormals)
  data(rtniIDs)
  rtniNormals<-get("rtniNormals")
  rtniIDs<-get("rtniIDs")
  if(sum(tnet%in%c("dpi","ref"))!=1)stop("'tnet' should be any of 'dpi' and 'ref'!")
  if(sum(idtype%in%c("entrez","probeid"))!=1)stop("'idtype' should be any of 'probeid' and 'entrez'!")
  if(eps>0)rtniNormals<-tni.dpi.filter(rtniNormals, eps=eps)
  if(idtype=="entrez"){
    tnetIDs=rtniIDs[,1:2]
  } else {
    tnetIDs=NULL
  }
  hits<-hits[!is.na(hits)]
  rtna <- tni2tna.preprocess(rtniNormals, hits=hits, verbose=FALSE)
  #compute regulons from both dpi and reference networks
  regs <- tna.mra(rtna, pValueCutoff=pValueCutoff, minRegulonSize=minRegulonSize, 
                  pAdjustMethod=pAdjustMethod, tnet=tnet, verbose=verbose)
  regs <- tna.get(regs,what="mra", ...=...)
  #remove regulons not listed in regref
  if(tnet=="dpi"){
    regref <- tna.mra(rtna, pValueCutoff=pValueCutoff, minRegulonSize=minRegulonSize, 
                      pAdjustMethod=pAdjustMethod, tnet="ref", verbose=verbose)
    regref <- tna.get(regref,what="mra", ...=...)
    idx<-rownames(regs)%in%rownames(regref)
    return(regs[idx,])
  } else {
    return(regs)
  }
}
##-----------------------------------------------------------------------------
Fletcher2013pipeline.mraTALL<-function(
  hits, minRegulonSize=20, idtype="probeid", pAdjustMethod="holm",
  tnet="dpi", eps=0, pValueCutoff=0.01, verbose=TRUE, ...){
  cat("Running analysis pipeline ... \n\n")
  data(rtniTALL)
  data(rtniIDs)
  rtniTALL<-get("rtniTALL")
  rtniIDs<-get("rtniIDs")
  if(sum(tnet%in%c("dpi","ref"))!=1)stop("'tnet' should be any of 'dpi' and 'ref'!")
  if(sum(idtype%in%c("entrez","probeid"))!=1)stop("'idtype' should be any of 'probeid' and 'entrez'!")
  if(eps>0)rtniTALL<-tni.dpi.filter(rtniTALL, eps=eps)
  if(idtype=="entrez"){
    tnetIDs=rtniIDs[,1:2]
  } else {
    tnetIDs=NULL
  }
  hits<-hits[!is.na(hits)]
  rtna <- tni2tna.preprocess(rtniTALL, hits=hits, verbose=FALSE)
  #compute regulons from both dpi and reference networks
  regs <- tna.mra(rtna, pValueCutoff=pValueCutoff, minRegulonSize=minRegulonSize, 
                  pAdjustMethod=pAdjustMethod, tnet=tnet, verbose=verbose)
  regs <- tna.get(regs,what="mra", ...=...)
  #remove regulons not listed in regref
  if(tnet=="dpi"){
    regref <- tna.mra(rtna, pValueCutoff=pValueCutoff, minRegulonSize=minRegulonSize, 
                      pAdjustMethod=pAdjustMethod, tnet="ref", verbose=verbose)
    regref <- tna.get(regref,what="mra", ...=...)
    idx<-rownames(regs)%in%rownames(regref)
    return(regs[idx,])
  } else {
    return(regs)
  }
}
##-----------------------------------------------------------------------------
Fletcher2013pipeline.masters<-function(){
  cat("Starting analysis in batch mode... \n\n")
  data(rtni1st)
  rtni1st<-get("rtni1st")
  #---cohort I
  mra1stExp1<-rownames(Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp1',idtype="entrez")))
  mra1stExp2<-rownames(Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp2',idtype="entrez")))
  mra1stExp3<-rownames(Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp3',idtype="entrez")))
  #---cohort II
  mra2ndExp1<-rownames(Fletcher2013pipeline.mra2nd(hits=get.hits(what='Exp1',idtype="entrez")))
  mra2ndExp2<-rownames(Fletcher2013pipeline.mra2nd(hits=get.hits(what='Exp2',idtype="entrez")))
  mra2ndExp3<-rownames(Fletcher2013pipeline.mra2nd(hits=get.hits(what='Exp3',idtype="entrez")))
  #---consensus
  Exp1mra<-intersect(mra1stExp1,mra2ndExp1)
  Exp2mra<-intersect(mra1stExp2,mra2ndExp2)
  Exp3mra<-intersect(mra1stExp3,mra2ndExp3)
  #---get masters
  masters<-unique(c(Exp1mra,Exp2mra,Exp3mra))
  idx<-rtni1st@annotation$ENTREZ%in%masters
  masters<-rtni1st@annotation[idx,]
  return(masters)
}
##-----------------------------------------------------------------------------
Fletcher2013pipeline.consensus<-function(){
  cat("Running analysis pipeline ... \n\n")
  data(miscellaneous)
  data(rtni1st)
  consensus<-get("consensus")
  rtni1st<-get("rtni1st")
  masters<-rtni1st@annotation[consensus$masters,"ENTREZ"]
  names(masters)<-rtni1st@annotation[consensus$masters,"SYMBOL"]
  rm(rtni1st)
  getcon<-function(masters,hits1,hits2,hits3,eps){
    #---cohort I
    mra1stExp1<-rownames(Fletcher2013pipeline.mra1st(hits=hits1, eps=eps))
    mra1stExp2<-rownames(Fletcher2013pipeline.mra1st(hits=hits2, eps=eps))
    mra1stExp3<-rownames(Fletcher2013pipeline.mra1st(hits=hits3, eps=eps))
    #---cohort II
    mra2ndExp1<-rownames(Fletcher2013pipeline.mra2nd(hits=hits1, eps=eps))
    mra2ndExp2<-rownames(Fletcher2013pipeline.mra2nd(hits=hits2, eps=eps))
    mra2ndExp3<-rownames(Fletcher2013pipeline.mra2nd(hits=hits3, eps=eps))
    #---Normals
    mraNormalsExp1<-rownames(Fletcher2013pipeline.mraNormals(hits=hits1, eps=eps))
    mraNormalsExp2<-rownames(Fletcher2013pipeline.mraNormals(hits=hits2, eps=eps))
    mraNormalsExp3<-rownames(Fletcher2013pipeline.mraNormals(hits=hits3, eps=eps)) 
    #---T-ALL
    mraTALLExp1<-rownames(Fletcher2013pipeline.mraTALL(hits=hits1,eps=eps))
    mraTALLExp2<-rownames(Fletcher2013pipeline.mraTALL(hits=hits2,eps=eps))
    mraTALLExp3<-rownames(Fletcher2013pipeline.mraTALL(hits=hits3,eps=eps))
    #---compute consensus
    a<-cbind(Exp1=masters%in%mra1stExp1,Exp2=masters%in%mra1stExp2,Exp3=masters%in%mra1stExp3)
    b<-cbind(Exp1=masters%in%mra2ndExp1,Exp2=masters%in%mra2ndExp2,Exp3=masters%in%mra2ndExp3)
    c<-cbind(Exp1=masters%in%mraNormalsExp1,Exp2=masters%in%mraNormalsExp2,Exp3=masters%in%mraNormalsExp3)
    d<-cbind(Exp1=masters%in%mraTALLExp1,Exp2=masters%in%mraTALLExp2,Exp3=masters%in%mraTALLExp3)
    tb<-cbind(a,b,c,d)
    rownames(tb)<-masters
    tb
  }
  #---get tables
  hits1=get.hits(what='Exp1',idtype="entrez")
  hits2=get.hits(what='Exp2',idtype="entrez")
  hits3=get.hits(what='Exp3',idtype="entrez")
  tb00<-getcon(masters,hits1,hits2,hits3,eps=0.00)
  tb01<-getcon(masters,hits1,hits2,hits3,eps=0.01)
  tb05<-getcon(masters,hits1,hits2,hits3,eps=0.05)
  #----------
  #---get consensustable  
  contab<-tb00*3
  contab[contab==0]<-as.numeric(tb01[contab==0])*2
  contab[contab==0]<-as.numeric(tb05[contab==0])
  #---reorder plot consensus
  tb<-contab
  rownames(tb)<-names(masters)
  idx<-sort.list(
    rowSums(tb[,1:6]>0)+
      rowSums(tb[,1:6]==3)*0.75+
      rowSums(tb[,1:6]==2)*0.25-
      rowSums(tb[,10:12]>0)*0.1,
    decreasing=FALSE)
  tb<-tb[rev(idx),]
  plot.summary.mra(tb)
  #return supermasters
  contab<-contab[rev(idx),]
  idx<-rowSums(contab[,1:6]>0)==6
  conmasters<-rownames(contab)[idx]
  conmasters<-masters[masters%in%conmasters]
  colnames(contab)<-paste(rep(c("C1","C2","NR","TA"),each=3),colnames(contab),sep=".")
  masters<-masters[match(rownames(contab),masters)]
  list(masters=masters, evidenceTable=contab, consensusmasters=conmasters)
}
Fletcher2013mra.consensus<-function(){
  data(miscellaneous)
  consensus<-get("consensus")
  masters<-consensus$masters
  #---reorder plot consensus
  tb<-consensus$evidenceTable
  rownames(tb)<-names(masters)
  idx<-sort.list(
    rowSums(tb[,1:6]>0)+
      rowSums(tb[,1:6]==3)*0.75+
      rowSums(tb[,1:6]==2)*0.25-
      rowSums(tb[,10:12]>0)*0.1,
    decreasing=FALSE)
  tb<-tb[rev(idx),]
  plot.summary.mra(tb)
}
##-----------------------------------------------------------------------------
Fletcher2013pipeline.siRNA<-function(what="siESR1"){
  cat("Running analysis pipeline ... \n\n")
  if(sum(what%in%c("siESR1","siSPDEF","siPTTG1"))!=1)stop("'what' should be any of 'siESR1','siSPDEF' and 'siPTTG1'!")
  data(rtni1st)
  data(siRNA)
  data(miscellaneous)
  rtni1st<-get("rtni1st")
  siRNA<-get("siRNA")
  consensus<-get("consensus")
  masters<-consensus$masters
  #--- get pheno e hits
  if(what=="siESR1"){
    hits<-siRNA$ESR1$ENTREZ
  } else if(what=="siSPDEF"){
    hits<-siRNA$SPDEF$ENTREZ
  } else if(what=="siPTTG1"){
    hits<-siRNA$PTTG1$ENTREZ
  }
  hits<-hits[!is.na(hits)]
  rtni1st@results$tn.dpi<-rtni1st@results$tn.dpi[,masters]
  rtni1st@results$tn.ref<-rtni1st@results$tn.ref[,masters]
  rtni1st@transcriptionFactors=masters
  rtna <- tni2tna.preprocess(rtni1st, hits=hits, verbose=FALSE)
  #compute regulons
  #ps. p.vaule is set to 1.0 to show all results in a table!!
  rtna <- tna.mra(rtna, pValueCutoff=1, minRegulonSize=20, pAdjustMethod="holm", verbose=TRUE)
  tna.get(rtna,what="mra")
}
##-----------------------------------------------------------------------------
get.hits<-function(what="Exp1", sig="FGFR",idtype="probeid"){
  if(sum(what%in%c("Exp1","Exp2","Exp3"))!=1)
    stop("'what' should be any of 'Exp1','Exp2' and 'Exp3'!")
  if(sum(sig%in%c("E2","FGFR","DT","Tet","random","uniqueE2","uniqueFGFR","intersectE2FGFR",'uniqueDT'))!=1)
    stop("'sig' should be any of 'E2','FGFR', 'DT', 'Tet', 
         'random', 'uniqueE2', 'uniqueFGFR', 'intersectE2FGFR' and 'uniqueDT'!")
  if(sig=="Tet" && what!="Exp3")
    stop("sig='Tet' only valid for what='Exp3'!")
  if(what=="Exp1"){
    sgt<-Fletcher2013pipeline.deg(what="Exp1",idtype=idtype)
    if(sig=="random"){
      hits=sgt$random
    } else if(sig=="DT"){
      hits=sgt$DT
    } else if(sig=="E2"){
      hits=sgt$E2
    } else if(sig=="FGFR"){
      hits=sgt$E2FGF10
    } else if(sig=="uniqueE2"){
      hits=setdiff(sgt$E2,sgt$E2FGF10)
    } else if(sig=="uniqueFGFR"){
      hits=setdiff(sgt$E2FGF10,sgt$E2)
    } else if(sig=="intersectE2FGFR"){
      hits=intersect(sgt$E2FGF10,sgt$E2)
    } else if(sig=="uniqueDT"){
      hits=setdiff(sgt$DT,unique(sgt$E2,sgt$E2FGF10))
    } else {
      stop("'get.hits' option not avilable!")
    }
  } else if(what=="Exp2"){
    sgt<-Fletcher2013pipeline.deg(what="Exp2",idtype=idtype)
    if(sig=="random"){
      hits=sgt$random
    } else if(sig=="DT"){
      hits=sgt$DT
    } else if(sig=="E2"){
      hits=sgt$E2      
    } else if(sig=="FGFR"){
      hits=sgt$E2AP20187      
    } else if(sig=="uniqueE2"){
      hits=setdiff(sgt$E2,sgt$E2AP20187)
    } else if(sig=="uniqueFGFR"){
      hits=setdiff(sgt$E2AP20187,sgt$E2)
    } else if(sig=="intersectE2FGFR"){
      hits=intersect(sgt$E2AP20187,sgt$E2)
    } else if(sig=="uniqueDT"){
      hits=setdiff(sgt$DT,unique(sgt$E2,sgt$E2AP20187))
    } else {
      stop("'get.hits' option not avilable!")
    }
  } else if(what=="Exp3"){
    sgt<-Fletcher2013pipeline.deg(what="Exp3",idtype=idtype)
    if(sig=="random"){
      hits=sgt$random
    } else if(sig=="DT"){
      hits=sgt$TetDT
    } else if(sig=="Tet"){
      hits=sgt$Tet      
    } else if(sig=="E2"){
      hits=sgt$TetE2      
    } else if(sig=="FGFR"){
      hits=sgt$TetE2FGF10
    } else if(sig=="uniqueE2"){
      hits=setdiff(sgt$TetE2,sgt$TetE2FGF10)
    } else if(sig=="uniqueFGFR"){
      hits=setdiff(sgt$TetE2FGF10,sgt$TetE2)
    } else if(sig=="intersectE2FGFR"){
      hits=intersect(sgt$TetE2FGF10,sgt$TetE2)
    } else if(sig=="uniqueDT"){
      hits=setdiff(sgt$TetDT,unique(sgt$TetE2,sgt$TetE2FGF10))
    } else {
      stop("'get.hits' option not avilable!")
    }
  }
  hits
}
