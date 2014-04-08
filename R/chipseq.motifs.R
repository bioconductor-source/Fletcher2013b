###############################################################################
# Motif analysis using bdsite maps derived from ChIP-seq experiments   				    
###############################################################################
Fletcher2013pipeline.chipseq<-function(what="SPDEF"){
  cat("Running analysis pipeline ... \n\n")
  opts<-c("ESR1","GATA3","FOXA1","SPDEF")
  if(sum(what%in%opts)!=1)stop(paste("'what' should be any of the options:",opts))
  #get regulons
  data(rtni1st,rtniIDs)
  data(miscellaneous)
  consensus<-get("consensus")
  rtni1st<-get("rtni1st")
  rtniIDs<-get("rtniIDs")
  ESR1bdsites<-get("ESR1bdsites")
  FOXA1bdsites<-get("FOXA1bdsites")
  GATA3bdsites<-get("GATA3bdsites")
  SPDEFbdsites<-get("SPDEFbdsites")
  chromlen<-get("chromlen")
  #get masters
  masters<-consensus$masters
  bin<-abs(rtni1st@results$tn.dpi[,masters[opts]])
  bin[bin!=0]=1
  #get TSS map
  tssmap<-getTSSmap(rtniIDs$ENTREZ)
  if(what=="ESR1"){
    #compute stats for ESR1 and plot
    gene.list<-getGeneDistBD(tssmap,ESR1bdsites$bdsites,chromlen)
    probe.list<-getProbeDistBD(rtniIDs, gene.list)
    bdpvals<-getstat(bin,probe.list,nper=1000)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="ESR1", ylim=c(0,0.030), file="ESR1chipESR1reg", adjust=1.5)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="GATA3", ylim=c(0,0.030), file="ESR1chipGATA3reg", adjust=1.5)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="FOXA1", ylim=c(0,0.030), file="ESR1chipFOXA1reg", adjust=1.5)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="SPDEF", ylim=c(0,0.030), file="ESR1chipSPDEFreg", adjust=1.5)
  } else if(what=="FOXA1"){
    #compute stats for FOXA1 and plot
    gene.list<-getGeneDistBD(tssmap,FOXA1bdsites$bdsites,chromlen)
    probe.list<-getProbeDistBD(rtniIDs, gene.list)
    bdpvals<-getstat(bin,probe.list,nper=1000)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="ESR1", ylim=c(0,0.030),file="FOXA1chipESR1reg", adjust=1.5)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="GATA3", ylim=c(0,0.030),file="FOXA1chipGATA3reg", adjust=1.5)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="FOXA1", ylim=c(0,0.030),file="FOXA1chipFOXA1reg", adjust=1.5)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="SPDEF", ylim=c(0,0.030),file="FOXA1chipSPDEFreg", adjust=1.5)
  } else if(what=="GATA3"){
    #compute stats for GATA3 and plot
    gene.list<-getGeneDistBD(tssmap,GATA3bdsites$bdsites,chromlen)
    probe.list<-getProbeDistBD(rtniIDs, gene.list)
    bdpvals<-getstat(bin,probe.list,nper=1000)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="ESR1", ylim=c(0,0.015),file="GATA3chipESR1reg")
    plotregs(probe.list, bdpvals, bin, masters, labmaster="GATA3", ylim=c(0,0.015),file="GATA3chipGATA3reg")
    plotregs(probe.list, bdpvals, bin, masters, labmaster="FOXA1", ylim=c(0,0.015),file="GATA3chipFOXA1reg")
    plotregs(probe.list, bdpvals, bin, masters, labmaster="SPDEF", ylim=c(0,0.015),file="GATA3chipSPDEFreg")
  } else if(what=="SPDEF"){  
    #compute stats for SPDEF and plot
    gene.list<-getGeneDistBD(tssmap,SPDEFbdsites$bdsites,chromlen)
    probe.list<-getProbeDistBD(rtniIDs, gene.list)
    bdpvals<-getstat(bin,probe.list,nper=1000)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="ESR1", ylim=c(0,0.015),file="SPDEFchipESR1reg", adjust=1.5)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="GATA3", ylim=c(0,0.015),file="SPDEFchipGATA3reg", adjust=1.5)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="FOXA1", ylim=c(0,0.015),file="SPDEFchipFOXA1reg", adjust=1.5)
    plotregs(probe.list, bdpvals, bin, masters, labmaster="SPDEF", ylim=c(0,0.015),file="SPDEFchipSPDEFreg", adjust=1.5)
  }
}
###############################################################################
# Motif analysis using pre-processed data derived FIMO + TRANSFAC PWMs				    
###############################################################################
Fletcher2013pipeline.motifs<-function(){
  cat("Running analysis pipeline ... \n\n")
  #get regulons
  data(rtni1st,rtniIDs)
  data(miscellaneous)
  consensus<-get("consensus")
  rtni1st<-get("rtni1st")
  rtniIDs<-get("rtniIDs")
  fimoESR1<-get("fimoESR1")
  fimoFOXA1<-get("fimoFOXA1")
  fimoGATA3<-get("fimoGATA3")
  masters<-consensus$masters
  bin<-abs(rtni1st@results$tn.dpi[,masters[c("ESR1","GATA3","FOXA1","SPDEF")]])
  bin[bin!=0]=1
  #------------ ESR1 regulon
  reg<-bin[,masters["ESR1"]]
  reg<-names(reg[reg==1])
  #compute stats for ESR1 and plot
  mts<-abs(fimoESR1$motifs[reg,"nearbd"])
  rdm<-abs(fimoESR1$rd[reg,"nearbd"])
  mts <-mts[!is.na(mts)]
  rdm <-rdm[!is.na(rdm)]
  ptest<-wilcox.test(mts, rdm, paired = FALSE, alternative = "less")
  plotallbd(fimoESR1, reg, labmaster="ESR1", ptest$p.value, ylim=c(0,0.020), file="ESR1fimo")
  #------------ GATA3 regulon
  reg<-bin[,masters["GATA3"]]
  reg<-names(reg[reg==1])
  #compute stats for GATA3 and plot
  mts<-abs(fimoGATA3$motifs[reg,"nearbd"])
  rdm<-abs(fimoGATA3$rd[reg,"nearbd"])
  mts <-mts[!is.na(mts)]
  rdm <-rdm[!is.na(rdm)]
  ptest<-wilcox.test(mts, rdm, paired = FALSE, alternative = "less")
  plotallbd(fimoGATA3, reg, labmaster="GATA3", ptest$p.value, ylim=c(0,0.020), file="GATA3fimo")
  #------------ FOXA1 regulon
  reg<-bin[,masters["FOXA1"]]
  reg<-names(reg[reg==1])
  #compute stats for FOXA1 and plot
  mts<-abs(fimoFOXA1$motifs[reg,"nearbd"])
  rdm<-abs(fimoFOXA1$rd[reg,"nearbd"])
  mts <-mts[!is.na(mts)]
  rdm <-rdm[!is.na(rdm)]
  ptest<-wilcox.test(mts, rdm, paired = FALSE, alternative = "less")
  plotallbd(fimoFOXA1, reg, labmaster="FOXA1", ptest$p.value, ylim=c(0,0.020), file="FOXA1fimo")
}


