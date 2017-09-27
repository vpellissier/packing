packing<-function(site.min, site.max, traits, mat, ncores=NULL) #here, site is the number/id of the cell
{
  
  convhulln(traits[site.min,], "FA")$vol->vol.min
  convhulln(traits[site.max,], "FA")$vol->vol.max
  
  setdiff(site.max, site.min)->esp.unique.max
  length(esp.unique.max)->nbsp.unique.max
  intersect(site.min, site.max)->esp.comm
  esp.rem<-NULL
  
  volume.depaup<-rich.depaup<-vector()
  
  l.a.max.depaup<-lapply(esp.unique.max, function(x) c(setdiff(esp.unique.max,x), esp.comm))
  
  for (i in seq(length(esp.unique.max)+length(esp.comm)))#length(esp.unique.max))
  {
    setdiff(esp.unique.max, esp.rem)->esp.unique.max
    
    l.a.max.depaup<-lapply(esp.unique.max, function(x) c(setdiff(esp.unique.max,x), esp.comm))
    
    if(!is.null(ncores))
    {
      force(l.a.max.depaup)
      sfExport("l.a.max.depaup", envir=environment())
      sfSapply(l.a.max.depaup, function(x) convhulln(traits[x,], options=c("FA", "QJ"))$vol)->fric.obs    
    }
    
    if(is.null(ncores))
    {
      sapply(l.a.max.depaup, function(x) convhulln(traits[x,], options=c("FA", "QJ"))$vol)->fric.obs
    }
    
    esp.unique.max[which.max(vol.max-fric.obs)]->esp.rem   #largest difference
    fric.obs[which.max(vol.max-fric.obs)]->volume.depaup[i]
    length(l.a.max.depaup[[which.max(vol.max-fric.obs)]])->rich.depaup[i]
    
    if(volume.depaup[i]<=vol.min) break
  }
  return(list(pourc.pack=(1-((i-1)/length(site.max)))*100, nb.esp.expansion=i-1, nb.esp.com=length(esp.comm), 
              nb.esp.min=length(site.min), nb.esp.max=length(site.max)))
}
#################################
# selection of 3 assemblages per biomes.
# 1/ most diverse assemblage at high NPP (A1= richest assemblage in the cells having high NPP, as defined in 9 classes)
# 2/ selection of an assemblage at medium NPP (A2= assemblage with 50% richness of A1, in NPP class 4)
# 3/ selection of an assemblage at low NPP (A2= assemblage with 50% richness of A2, in NPP class 1)


pack.biomes<-function(m, niter=10, ncores=NULL, tab, traits)
{
  packs<-data.frame(matrix(ncol=20, nrow=0))
  names(packs)<-c("RxB", "Sample", "NPP low", "NPP medium", "NPP high",
                  paste(c("pourc.pack", "nb.exp.expansion", "nb.esp.com", "nb.esp.min", "nb.esp.max"),
                        rep(c("lh", "mh", "lm"), each=5), sep="."))
  
  for (i in seq(niter))
  {
    tab[tab$realms_x_biomes==m,]->rs.tab.region
    if(m=="R1.B2") 
      rs.tab.region<-rs.tab.region[rs.tab.region$nppt<400,]
    
    rs.tab.region$cell<-as.character(rs.tab.region$cell)
    cooc[is.element(rownames(cooc), rs.tab.region$cell),]->cooc_region
    rs.tab.region$class.npp<-Hmisc::cut2(rs.tab.region$nppt, g=9)
    rs.tab.region$class.npp<-cut(rs.tab.region$nppt, breaks=9, labels=F)
    rs.tab.region$labels.class.npp<-cut(rs.tab.region$nppt, breaks=9)
    
    if(length(unique(rs.tab.region$class.npp))!=9) 
      return(packs)
    
    rs.tab.region.high<-rs.tab.region[rs.tab.region$class.npp==9,]
    rs.tab.region.medium<-rs.tab.region[rs.tab.region$class.npp==5,]
    rs.tab.region.low<-rs.tab.region[rs.tab.region$class.npp==1,]
    
    rs.tab.region.high<-rs.tab.region.high[rs.tab.region.high$nbsp>max(rs.tab.region.high$nbsp)*.90,]    #high cell sampled in cells within 95% of max nbsp
    
    cell.high.npp<-sample(rs.tab.region.high$cell,1)
    nbsp.high.npp<-rs.tab.region.high$nbsp[rs.tab.region.high$cell==cell.high.npp]
    
    cell.medium<-rs.tab.region.medium[which.min(abs(rs.tab.region.medium$nbsp-nbsp.high.npp/2)), "cell"] #cell closest to x% of the high npp cell
    nbsp.medium<-rs.tab.region.medium[rs.tab.region.medium$cell==cell.medium, "nbsp"]
    
    cell.medium.npp<-sample(rs.tab.region.medium[rs.tab.region.medium$nbsp>nbsp.medium*0.95 & rs.tab.region.medium$nbsp<nbsp.medium*1.05, "cell"],1)  #sampling one within all cells having ±5% species than the selected cell
    nbsp.medium.npp<-rs.tab.region.medium[rs.tab.region.medium$cell==cell.medium.npp, "nbsp"]
    
    rm(list=c("cell.medium", "nbsp.medium"))
    
    cell.low<-rs.tab.region.low[which.min(abs(rs.tab.region.low$nbsp-nbsp.medium.npp/2)), "cell"] #cell closest to x% of the high npp cell
    nbsp.low<-rs.tab.region.low[rs.tab.region.low$cell==cell.low, "nbsp"]
    
    cell.low.npp<-sample(rs.tab.region.low[rs.tab.region.low$nbsp>nbsp.low*0.95 & rs.tab.region.low$nbsp<nbsp.low*1.05, "cell"],1)  #sampling one within all cells having ±5% species than the selected cell
    nbsp.low.npp<-rs.tab.region.low[rs.tab.region.low$cell==cell.low.npp, "nbsp"]
    
    rm(list=c("cell.low", "nbsp.low"))
    
    if(!is.null(ncores)) 
    {
      sfInit(parallel = TRUE, cpus=40)
      sfExport("traits")
      sfLibrary(geometry)
    }
    
    cell.low.npp<-cooc_region[cell.low.npp,]
    cell.medium.npp<-cooc_region[cell.medium.npp,]
    cell.high.npp<-cooc_region[cell.high.npp,]
    
    cell.low.npp<-names(cell.low.npp[cell.low.npp>0])
    cell.medium.npp<-names(cell.medium.npp[cell.medium.npp>0])
    cell.high.npp<-names(cell.high.npp[cell.high.npp>0])
    
    packing(cell.low.npp, cell.high.npp, traits=traits, mat=cooc_region, ncores=ncores)->pack.LH
    packing(cell.medium.npp, cell.high.npp, traits=traits, mat=cooc_region, ncores=ncores)->pack.MH
    packing(cell.low.npp, cell.medium.npp, traits=traits, mat=cooc_region, ncores=ncores)->pack.LM
    
    if(!is.null(ncores)) sfStop()
    
    as.character(unique(rs.tab.region$labels.class.npp[rs.tab.region$class.npp==1]))->npp.l
    as.character(unique(rs.tab.region$labels.class.npp[rs.tab.region$class.npp==4]))->npp.m
    as.character(unique(rs.tab.region$labels.class.npp[rs.tab.region$class.npp==9]))->npp.h
    
    data.frame(t(c(m, i, npp.l, npp.m, npp.h, unlist(pack.LH), unlist(pack.MH), unlist(pack.LM))))->pack.region
    
    names(packs)->names(pack.region)
    
    rbind(packs, pack.region)->packs
    print(i)
  }
  print(paste(m, "Try n°", i))
  return(packs)
}