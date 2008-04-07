
# a function to color bait nodes yellow and prey nodes white

# colors reciprocated edges red and put arrows on both ends

# colors unreciprocated edges connecting baits blue and puts an arrow pointing
# toward the bait that was detected but did not detect

#colors unreciprocated edges from baits to prey gray and puts an arrow
# pointing toward the prey

plotComplex =
function(complexMembers,g,VBs,VPs,geneName=FALSE,baitColor="yellow",preyColor="white",recipLineColor="red",unrecipBBLineColor="blue",unrecipBPLineColor="gray",y="neato"){

  ##tc says..most of the graphs that i create are intactGraphs which extends graphs
  #stopifnot(class(g)=="graphNEL")
  
  
  sg=subGraph(complexMembers,g)
  
  nAttrs = list()
  eAttrs = list()
  
  baits = intersect(complexMembers,VBs)
  
  nAttrs$fillcolor = rep(preyColor,length(complexMembers))
  names(nAttrs$fillcolor) = complexMembers
  nAttrs$fillcolor[baits] = baitColor
  
  if(geneName){
	nodeLabels = unlist(mget(complexMembers,env=YEASTGENENAME))
	nAttrs$label = nodeLabels
	}



  
  eNames = edgeNames(sg)
  
  eAttrs$arrowhead = rep("none",length(eNames))
  names(eAttrs$arrowhead) = eNames
  eAttrs$arrowtail = rep("none",length(eNames))
  names(eAttrs$arrowtail) = eNames	
  eAttrs$color = rep(unrecipBPLineColor,length(eNames))
  names(eAttrs$color) = eNames
  
  for (k in eNames){
    
    twonodes = unlist(strsplit(k,"~"))
    numBs = sum(twonodes %in% VBs)
    
    if(numBs==1){
      bnode = which(twonodes %in% VBs)
      if(bnode==1){
        eAttrs$arrowtail[k] = "open"
      } else eAttrs$arrowhead[k] = "open"
     }
    
    if(numBs==2){
      ssg=subGraph(twonodes,sg)
      
      if(numEdges(ssg)==2){
        eAttrs$arrowhead[k]="open"
        eAttrs$arrowtail[k]="open"
        eAttrs$color[k]=recipLineColor
      }
      
      if(numEdges(ssg)==1){
        eAttrs$color[k]=unrecipBBLineColor
        outnode = names(which(degree(ssg)$outDegree==1))
        intwo = which(twonodes==outnode)
        if(intwo==1){
          eAttrs$arrowtail[k]="open"
        } else eAttrs$arrowhead[k]="open"	
      }
    }
  }
  
  plot(sg,nodeAttrs=nAttrs,edgeAttrs=eAttrs,y=y)
  
}
