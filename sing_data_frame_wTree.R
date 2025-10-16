singleton.to.data.frame.withoutTree <- function(spp, lat, long){
  names(lat) <- spp
  names(long) <- spp 
  locais <- cbind(lat, long)
  
  singletons <- NULL
  sppsORD <- sort(rownames(locais), decreasing = F)
  for(i in 1:length(unique(rownames(locais)))){
    singletons[i] <- length(grep(pattern = unique(sppsORD)[i], x = sppsORD))
  }
  
  if(all(singletons > 1) == TRUE){
    spp <- spp
    data_df <- data.frame(spp = spp, long = long, lat = lat)
  } else if(any(singletons == 1) == TRUE){
    singletons <- which(singletons == 1)
    um_ponto <- unique(sppsORD)[singletons]
    locaisNew <- locais
    for(i in 1:length(um_ponto)){
      linha <- which(rownames(locaisNew) == um_ponto[i])
      locaisNew <- locaisNew[-linha,]
    }
    
    # data frame
    lat <- locaisNew[, 1]
    long <- locaisNew[, 2]
    spp <- rownames(locaisNew)
    
    data_df <- data.frame(spp = spp, long = long, lat = lat)
  }
}