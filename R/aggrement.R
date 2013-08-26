###############################################################################
#Plot aggrement among signatures    			    
###############################################################################
Fletcher2013mra.agreement.cohort1<-function(){
  cat("Starting analysis in batch mode... \n\n")
  Exp1<-Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp1',idtype="entrez"), tnet="ref", minRegulonSize=20, ntop=-1, order=FALSE, reportNames=FALSE)
  Exp2<-Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp2',idtype="entrez"), tnet="ref", minRegulonSize=3,ntop=-1, order=FALSE, reportNames=FALSE)
  Exp3<-Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp3',idtype="entrez"), tnet="ref", minRegulonSize=3, ntop=-1, order=FALSE, reportNames=FALSE)
  xrand<-Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp1',idtype="entrez",sig="random"), tnet="ref", minRegulonSize=3, ntop=-1, order=FALSE, reportNames=FALSE)
  sublist<-getmatch(lts=list(Exp1, Exp2, Exp3, xrand),glabs=c("E2.FGF10","E2.AP20187","Tet.E2.FGF10","Random"))
  pdf('agree.experiments.cohortI.pdf',6.5,6.8)
  plotcon(sublist, tlt="", glabs=c("Exp1","Exp2","Exp3","Random"),
          cols=c("red","blue","green","grey50"), top=50,cexp=0.35)
  dev.off()
  cat("-File 'agree.experiments.cohortI.pdf' generated!\n\n")
}
##-----------------------------------------------------------------------------
Fletcher2013mra.agreement.exp1<-function(){
  cat("Starting analysis in batch mode... \n\n")
  Exp1.1st<-Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp1',idtype="entrez"), tnet="ref", minRegulonSize=20, ntop=-1, order=FALSE,reportNames=FALSE)
  Exp1.2nd<-Fletcher2013pipeline.mra2nd(hits=get.hits(what='Exp1',idtype="entrez"), tnet="ref", minRegulonSize=3,ntop=-1, order=FALSE,reportNames=FALSE)
  Exp1.Norm<-Fletcher2013pipeline.mraNormals(hits=get.hits(what='Exp1',idtype="entrez"), tnet="ref", minRegulonSize=3, ntop=-1, order=FALSE,reportNames=FALSE)
  Exp1.TALL<-Fletcher2013pipeline.mraTALL(hits=get.hits(what='Exp1',idtype="entrez"), tnet="ref", minRegulonSize=3, ntop=-1, order=FALSE,reportNames=FALSE)
  sublist<-getmatch(lts=list(Exp1.1st, Exp1.2nd, Exp1.Norm, Exp1.TALL),glabs=c("Exp1.1st","Exp1.2nd","Exp1.Norm","Exp1.TALL"))
  pdf('agree.cohorts.exp1.pdf',6.5,6.8)
  plotcon(sublist, tlt="", glabs=c("Cohort I","Cohort II","Normals","T-ALL"),
          cols=colorRampPalette(c("darkred","grey50"))(4), top=50, cexp=0.35)
  dev.off()
  cat("-File 'agree.cohorts.exp1.pdf' generated!\n\n")
}
##-----------------------------------------------------------------------------
Fletcher2013mra.agreement.exp2<-function(){
  cat("Starting analysis in batch mode... \n\n")
  Exp2.1st<-Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp2',idtype="entrez"), tnet="ref", minRegulonSize=20, ntop=-1, order=FALSE,reportNames=FALSE)
  Exp2.2nd<-Fletcher2013pipeline.mra2nd(hits=get.hits(what='Exp2',idtype="entrez"), tnet="ref", minRegulonSize=3,ntop=-1, order=FALSE,reportNames=FALSE)
  Exp2.Norm<-Fletcher2013pipeline.mraNormals(hits=get.hits(what='Exp2',idtype="entrez"), tnet="ref", minRegulonSize=3, ntop=-1, order=FALSE,reportNames=FALSE)
  Exp2.TALL<-Fletcher2013pipeline.mraTALL(hits=get.hits(what='Exp2',idtype="entrez"), tnet="ref", minRegulonSize=3, ntop=-1, order=FALSE,reportNames=FALSE)
  sublist<-getmatch(lts=list(Exp2.1st, Exp2.2nd, Exp2.Norm, Exp2.TALL),glabs=c("Exp2.1st","Exp2.2nd","Exp2.Norm","Exp2.TALL"))
  pdf('agree.cohorts.exp2.pdf',6.5,6.8)
  plotcon(sublist, tlt="", glabs=c("Cohort I","Cohort II","Normals","T-ALL"),
          cols=colorRampPalette(c("darkblue","grey50"))(4), top=50, cexp=0.35)
  dev.off()
  cat("-File 'agree.cohorts.exp2.pdf' generated!\n\n")
}
##-----------------------------------------------------------------------------
Fletcher2013mra.agreement.exp3<-function(){
  cat("Starting analysis in batch mode... \n\n")
  Exp3.1st<-Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp3',idtype="entrez"), tnet="ref", minRegulonSize=20, ntop=-1, order=FALSE,reportNames=FALSE)
  Exp3.2nd<-Fletcher2013pipeline.mra2nd(hits=get.hits(what='Exp3',idtype="entrez"), tnet="ref", minRegulonSize=3,ntop=-1, order=FALSE,reportNames=FALSE)
  Exp3.Norm<-Fletcher2013pipeline.mraNormals(hits=get.hits(what='Exp3',idtype="entrez"), tnet="ref", minRegulonSize=3, ntop=-1, order=FALSE,reportNames=FALSE)
  Exp3.TALL<-Fletcher2013pipeline.mraTALL(hits=get.hits(what='Exp3',idtype="entrez"), tnet="ref", minRegulonSize=3, ntop=-1, order=FALSE,reportNames=FALSE)
  sublist<-getmatch(lts=list(Exp3.1st, Exp3.2nd, Exp3.Norm, Exp3.TALL),glabs=c("Exp3.1st","Exp3.2nd","Exp3.Norm","Exp3.TALL"))
  pdf('agree.cohorts.exp3.pdf',6.5,6.8)
  plotcon(sublist, tlt="", glabs=c("Cohort I","Cohort II","Normals","T-ALL"),
          cols=colorRampPalette(c("darkgreen","grey50"))(4), top=50, cexp=0.35)
  dev.off()
  cat("-File 'agree.cohorts.exp3.pdf' generated!\n\n")
}
