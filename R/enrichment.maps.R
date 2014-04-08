###############################################################################
#Return the enrichmen map in an igraph object      		    
###############################################################################
Fletcher2013pipeline.enrichmap<-function(){
  cat("Starting analysis in batch mode... \n\n")
  #get FGFR2 signatures in rtni1st
  mra1stX1<-Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp1',idtype="entrez"), tnet="ref", reportNames=FALSE)
  mra1stX2<-Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp2',idtype="entrez"), tnet="ref", reportNames=FALSE)
  mra1stX3<-Fletcher2013pipeline.mra1st(hits=get.hits(what='Exp3',idtype="entrez"), tnet="ref", reportNames=FALSE)
  #get results
  list.mra1st<-list()
  list.mra1st$exp1$FGF10<-rownames(mra1stX1)
  list.mra1st$exp2$FGF10<-rownames(mra1stX2)
  list.mra1st$exp3$FGF10<-rownames(mra1stX3)
  #------------
  #Call RedeR
  rdp<-RedPort()
  calld(rdp,maxlag=5000)
  resetd(rdp)
  cat("\nComputing network... \n\n")
  #------------
  data(rtni1st)
  data(miscellaneous)
  rtni1st<-get("rtni1st")
  consensus<-get("consensus")
  masters<-consensus$masters
  #------------
  list.mra1st$exp1$FGF10<-rownames(rtni1st@annotation)[rtni1st@annotation$ENTREZ%in%list.mra1st$exp1$FGF10]
  list.mra1st$exp2$FGF10<-rownames(rtni1st@annotation)[rtni1st@annotation$ENTREZ%in%list.mra1st$exp2$FGF10]
  list.mra1st$exp3$FGF10<-rownames(rtni1st@annotation)[rtni1st@annotation$ENTREZ%in%list.mra1st$exp3$FGF10]
  #------------
  #compute jaccard coefficient (filter JC < 0.4)
  #obs1. full-overlap is shown in heatmap figures!
  tn.ref<-abs(rtni1st@results$tn.ref)
  tn.ref<-tn.ref[,colSums(tn.ref!=0)>=20]
  adjmt<-getJC(tn.ref)
  adjmt[adjmt<0.40]=0
  #get regulon's degree
  bin<-tn.ref
  bin[bin>0]=1
  degree<-colSums(bin)
  ##add symbols
  tfs<-tni.get(rtni1st, what = "tfs")
  idx<-tfs%in%colnames(adjmt)
  masterids<-data.frame(probeid=tfs[idx],symbol=names(tfs[idx]))
  rownames(masterids)<-tfs[idx]
  #add general hits
  masterids<-data.frame(masterids,Masters=as.numeric(rownames(masterids)%in%masters))
  masterids<-data.frame(masterids,Genral1=as.numeric(rownames(masterids)%in%list.mra1st$exp1$FGF10))
  masterids<-data.frame(masterids,Genral2=as.numeric(rownames(masterids)%in%list.mra1st$exp2$FGF10))
  masterids<-data.frame(masterids,Genral3=as.numeric(rownames(masterids)%in%list.mra1st$exp3$FGF10))
  hits<-rowSums(masterids[,c("Genral1","Genral2","Genral3")])
  masterids<-data.frame(masterids, HitsGen = hits)
  #add degree
  masterids<-data.frame(masterids, degree=degree[rownames(masterids)]) 
  #make graph
  g<-graph.adjacency(adjmt, diag=FALSE, mode="undirected", weighted=TRUE)
  g<-att.mapv(g=g,dat=masterids,refcol=1)
  g<-att.setv(g=g, from="symbol", to='nodeAlias')
  g<-att.sete(g=g, from="weight", to='edgeWidth', nquant=10, xlim=c(5,100,1))
  g<-att.setv(g=g, from="degree", to='nodeSize', xlim=c(50,210,1), nquant=10, roundleg=1)
  V(g)$nodeFontColor="black"
  E(g)$linkType="notnested"
  #set attributes
  mypalette<-colorRampPalette(brewer.pal(9,"Purples"), bias=1)(100)
  E(g)$edgeColor<-mypalette[30]
  V(g)$nodeLineColor<-mypalette[30]
  cols<-c(mypalette[10],colorRampPalette(c("white","orange","darkred"))(14)[c(3,5,9)])
  g<-att.setv(g=g, from="HitsGen", to='nodeColor',cols=cols)
  #cols<-c("#e4e4ff",colorRampPalette(c("white","orange","darkred"))(14)[c(3,5,9)])
  #g<-att.setv(g=g, from="HitsGen", to='nodeColor',cols=cols)
  V(g)$nodeFontSize<-80
  V(g)$nodeFontSize[V(g)$HitsGen==1]=80
  V(g)$nodeFontSize[V(g)$HitsGen>=2]=100
  #--add graph
  addGraph(rdp, g, layout=NULL, zoom=10)
  #--add legends
  addLegend.color(rdp,g, position="topright", title="FGFR2 enrichment", dxtitle=15, bend=0.6, ftsize=7,dxborder=10)
  addLegend.size(rdp,g,position="topleft", vertical=TRUE, title="Regulon size", dxtitle=36, ftsize=7, dxborder=10, dyborder=10)
  addLegend.size(rdp,g,position="topleft", type="edge", title="Jaccard coefficient", dxtitle=40, ftsize=7, dyborder=300, 
                 dxborder=10, intersp=40, edgelen=150)
  #--open layout dialog box
  relax(rdp,300,50,100,ps=TRUE)
  #return(g)
}

