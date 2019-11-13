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
