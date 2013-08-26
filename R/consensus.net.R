###############################################################################
#Return the enrichmen map in an igraph object        	    
###############################################################################
Fletcher2013pipeline.consensusnet<-function(){
  cat("Running analysis pipeline ... \n\n")
  data(rtni1st,rtniIDs)
  rtni1st<-get("rtni1st")
  rtniIDs<-get("rtniIDs")
  data(miscellaneous)
  consensus<-get("consensus")
  #--- get tfs
  tfs<-rtni1st@transcriptionFactors
  #--- get masters
  consensusmasters<-consensus$consensusmasters
  #--- get regulons
  rtni<-tni.dpi.filter(rtni1st, eps=0.01)
  adjmt<-getmat(consensusmasters,rtni@results$tn.dpi)
  #--- get signatures
  X1<- Fletcher2013pipeline.deg(what="Exp1")
  X2<- Fletcher2013pipeline.deg(what="Exp2")
  X3<- Fletcher2013pipeline.deg(what="Exp3")
  allSignatures<-unique(c(X1$E2FGF10,X2$E2AP20187,X3$TetE2FGF10))
  #------------
  #Call RedeR
  rdp<-RedPort()
  calld(rdp,maxlag=5000)
  resetd(rdp)
  cat("\nComputing network... \n\n")
  #------------
  #--- get ids, etc.!
  idx<-rownames(rtniIDs)%in%rownames(adjmt)
  masterids<-rtniIDs[idx,]
  masterids<-data.frame(masterids,TFs=as.numeric(rownames(masterids)%in%tfs))
  masterids<-data.frame(masterids,Masters=as.numeric(rownames(masterids)%in%consensusmasters))
  masterids<-data.frame(masterids,Exp1=as.numeric(rownames(masterids)%in%X1$E2FGF10))
  masterids<-data.frame(masterids,Exp2=as.numeric(rownames(masterids)%in%X2$E2AP20187))
  masterids<-data.frame(masterids,Exp3=as.numeric(rownames(masterids)%in%X3$TetE2FGF10))
  masterids<-data.frame(masterids,allSignatures=as.numeric(rownames(masterids)%in%allSignatures))
  sighits<-rowSums(masterids[,c("Exp1","Exp2","Exp3")])
  masterids<-data.frame(masterids, sighits=sighits)
  temp<-rep(0,nrow(masterids))
  temp[(masterids$Exp1)==1]=1
  temp[masterids$Masters==1]=2
  masterids<-data.frame(masterids,SigExp1=temp)
  temp<-rep(0,nrow(masterids))
  temp[(masterids$Exp2)==1]=1
  temp[masterids$Masters==1]=2
  masterids<-data.frame(masterids,SigExp2=temp)
  temp<-rep(0,nrow(masterids))
  temp[(masterids$Exp3)==1]=1
  temp[masterids$Masters==1]=2
  masterids<-data.frame(masterids,SigExp3=temp)
  temp<-rep(0,nrow(masterids))
  temp[(masterids$sighits)==3]=1
  temp[masterids$Masters==1]=2
  masterids<-data.frame(masterids,SigAllMas_012=temp)
  #--- get graph
  g<-graph.adjacency(adjmt, diag=FALSE, mode="undirected", weighted=TRUE)
  #--- load g attributes
  g<-att.mapv(g=g,dat=masterids,refcol=1)
  g<-att.setv(g=g, from="SYMBOL", to='nodeAlias')
  edl<-as.data.frame(get.edgelist(g),stringsAsFactors=FALSE)
  direction<-as.numeric(edl$V1%in%tfs)+as.numeric(edl$V2%in%tfs)
  direction[direction==2]=0
  E(g)$arrowDirection=direction
  E(g)$linkType="notnested"
  #--- set g attributes
  cols<-c(colorRampPalette(c("#e4e4ff","orchid","darkorchid4"))(10)[c(1,4,6)],"darkorchid4")
  g<-att.setv(g=g, from="sighits", to='nodeColor',cols=cols)
  g<-att.setv(g=g, from="sighits", to='nodeLineColor',cols=cols)
  g<-att.setv(g=g, from="Masters", to='nodeFontColor',cols=c("black","black"))
  g<-att.setv(g=g, from="TFs", to='nodeShape')
  g<-att.setv(g=g, from="Masters", to='nodeSize', xlim=c(55,120,1))
  V(g)$nodeFontSize<-V(g)$SigAllMas_012^1.3*50+1
  E(g)$edgeWidth<-3
  E(g)$edgeWidth<-3
  E(g)$edgeColor<-"grey80"
  E(g)$arrowLength=80
  E(g)$arrowAngle=15
  #--add graph
  addGraph(rdp, g, layout=NULL, zoom=10)
  #--add legends
  addLegend.color(rdp,g,title="FGFR2 enrichment (hits)", dxtitle=12, bend=0.6, ftsize=7,dxborder=15)
  addLegend.shape(rdp,g,title="", dxtitle=40, labvec=c("nTFs","TFs"), dyborder=70, dxborder=28, ftsize=9)
  relax(rdp,p1=90,ps=TRUE)
  #--- return graph
  #selectNodes(rdp,V(g)$name[V(g)$Masters==1])
  #return(g)
}
