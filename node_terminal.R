### Nova função para estimar MST de cada espécie de um dado nó (ou nós), tanto internos como terminais
# arguments:
# -coordin => coordinates
# -tree => phylogenetic tree
# -shape_file => polygon of some place already with datum
# -resol => resolution of quadrats (number of quadrats)
# -seeres => to see the resolution on a map
# -sobrepo => plotting onto a map already produced
# -mintreeall => just one minimum spanning tree (mst) of all species 
# -caption => legend
# -nodes => node of a tree with descendents
# -taxon => species in order to estimate msn graphs
# -poly => minimum convex polygon
# -transp => transparency degree of polygon
# -seephylog => if you want to see the phylogeny with nodes onto the map
terminal_node <- function(coordin, tree = NULL, shape_file, resol, seeres = FALSE, sobrepo = FALSE, mintreeall = FALSE, caption = TRUE,
                          nodes = NULL, taxon = NULL, pol = FALSE, transp = 0.5, seephylog = FALSE,
                          xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL){
  
  print('1) loading functions...')
  library(geiger); library(devtools)
  library(phytools)
  library(maptools); library(letsR); library(fossil)
  library(viridis); library(mapdata); library(vegan); library(spdep)
  library(dismo); library(spatstat); library(terra)
  
  if(length(taxon) < 2){
    warning('It is not possible with only one species to produce generalized tracks by the PAE-PCE analysis!')
  }
  
  if(sobrepo == FALSE){
    layout(1, 1, 1)
    plot(shape_file, type = 'n', xlim = c(xmin, xmax), ylim = c(ymin, ymax))
    # dev.off()
  }

  print('2) preparing the coordinates...')
  coordin <- matrix(as.matrix(coordin[, c(2, 3)]), nrow(coordin),
             2, dimnames = list(coordin[,1], colnames(coordin)[c(2, 3)]))
  
  cols1 <- setNames(viridis(n = length(unique(rownames(coordin)))),
                    unique(rownames(coordin)))
  
  if(seephylog == TRUE && sobrepo == TRUE){
    stop('This is not possible! Please try to use the each argument separately.')
  }


  #- R graphs:
  par(pty = 's')
  if(seephylog == TRUE){
    layout(1:2, 1, 2)
  }
  if(seeres == TRUE){
    layout(matrix(1, nr = 1, nc = 1, byrow = F))
  }
  if(seephylog == TRUE && seeres == TRUE){
    layout(matrix(c(1, 2, 3), nrow = 3, ncol = 1, byrow = F))
  }
  
  if(!is.null(tree)){
    print('3) preparing the tree...')
    # phylogenetic trees and nodes:
    if(is.null(tree$edge.length) == TRUE){
      print('...creating the branch lengths of a tree equal to one...')
      tree <- compute.brlen(tree, 1)
    }
    
    print('...checking names from dataset and names in the species tree.')
    if(name.check(tree, coordin)[1] != 'OK'){
      chk <- name.check(tree, coordin)
      tree <- drop.tip(tree, chk$tree_not_data)
    }
  }
  
  # plotting the tree and coordinates:
  if(is.na(proj4string(as(shape_file, 'Spatial'))) == TRUE){
    print('Attention: your shapefile has not been associated to any datum and this routine is going
            to associate it to WGS84 datum!')
    # proj4string(shape_file) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
    # projection(as(shape_file, 'Spatial')) <- CRS("+proj=longlat +datum=WGS84")
    crs(shape_file) <- "+proj=longlat +datum=WGS84"
  } else if(proj4string(as(shape_file, 'Spatial')) != "+proj=longlat +datum=WGS84"){
    # proj4string(shape_file) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
    # projection(as(shape_file, 'Spatial')) <- CRS("+proj=longlat +datum=WGS84")
    crs(shape_file) <- "+proj=longlat +datum=WGS84"
  }
  
  if(seephylog == TRUE){
    obj2 <- phylo.to.map(tree, coordin[,2:1], database = shape_file, plot = F, rotate = F)
    plot(obj2, colors = cols1, ftype = 'i', fsize = 0.6, cex.points = c(0.7, 1.2), pts = F,
         direction = "rightwards")
    labelnodes(1:(Ntip(tree) + tree$Nnode), 1:(Ntip(tree) + tree$Nnode), interactive = F, cex = .6)
    mtext('Phylogeny on the map', side = 3, line = 1)
  }
  
  if(!is.null(tree)){
    print('4) dividing polygon into grid and producing a raster...')
  } else {
    print('3) dividing polygon into grid and producing a raster...')
  }
  
  if(!is.null(tree)){
    arv <- reorder(tree, "postorder") # reordering the levels
    e1 <- arv$edge[, 1] # internal nodes
    e2 <- arv$edge[, 2] # terminal nodes and root
    EL <- arv$edge.length # branch lengths
    
    tabelao <- mrca(phy = tree, full = F)
  }
  

  grid <- raster(extent(as(shape_file, 'Spatial')), resolution = resol,
                 crs = CRS("+proj=longlat +datum=WGS84"))
  grid <- raster::extend(grid, c(1, 1))
  gridPolygon <- rasterToPolygons(grid)
  # suppressWarnings(proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84")) # datum WGS84
  #proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  # projection(gridPolygon) <- CRS("+proj=longlat +datum=WGS84")
  
  
  # clipping the intersected cells:
  # suppressWarnings(cropped_map <- raster::intersect(gridPolygon, shape_file))
  cropped_map <- raster::intersect(gridPolygon, as(shape_file, 'Spatial'))
  if (seeres == TRUE){
    plot(cropped_map, xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = T)
    mask.raster <- raster(extent(as(shape_file, 'Spatial')), resolution = resol,
                          crs = CRS("+proj=longlat +datum=WGS84"))
    r <- rasterize(as(shape_file, 'Spatial'), mask.raster, fun = 'first')
    proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
    r <- merge(r, mask.raster)
    ncellR <- ncell(r)
    for(j in 1:ncellR){
      if(is.na(r[j]) == FALSE){
        values(r)[j] <- 1
      }
    }
    map.r <- as.data.frame(rasterToPoints(r))
    pontosRaster <- rasterize(cbind(map.r$x, map.r$y), r, field = 1)
    # plot(pontosRaster, add = T)
    map.r$gridNumber <- which(pontosRaster@data@values == 1)
    text(map.r[,c(1, 2)], labels = map.r$gridNumber, cex = 0.8,
         col = adjustcolor(col = 'red', alpha = 1), font = 2)
  }
  
  # producing a raster of the shapefile
  mask.raster <- raster(extent(as(shape_file, 'Spatial')), resolution = resol,
                        crs = CRS("+proj=longlat +datum=WGS84"))
  suppressWarnings(r <- rasterize(as(shape_file, 'Spatial'), mask.raster))
  proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  # mask.raster[is.na(mask.raster)] <- 0
  r <- merge(r, mask.raster)
  ncellras <- ncell(r)
  
  dir.create('out/') # temporary folder
  
  if(sobrepo == FALSE){
    if(seephylog == TRUE){
      if(is.null(nodes)){
        plot(shape_file, axes = TRUE, sub = paste0(c('resolution: ', resol[1],
                    ' x ', resol[2]), collapse = ''), xlim = c(xmin, xmax), ylim = c(ymin, ymax))
        abline(h = 0, col = 'red', lty = 2) # equator
        mtext('Minimum spanning trees of the terminal nodes', side = 3, line = 1)
      } else if(!is.null(nodes)){
        if(length(nodes) == 1){
          plot(shape_file, axes = TRUE, sub = paste0(c('resolution: ', resol[1], ' x ', resol[2]),
                  collapse = ''), cex.main = 0.7, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
          abline(h = 0, col = 'red', lty = 2) # equator
          mtext(paste0(c('Mapping a minimum spanning tree of the internal node:',
                       nodes), collapse = ' '), side = 3, line = 1)
        } else if(length(nodes) > 1){
          plot(shape_file, axes = TRUE, sub = paste0(c('resolution: ',
                resol[1], ' x ', resol[2]), collapse = ''), cex.main = 0.7,
               xlim = c(xmin, xmax), ylim = c(ymin, ymax))
          abline(h = 0, col = 'red', lty = 2) # equator
          mtext(paste0(c('Mapping a minimum spanning tree of the internal nodes:',
                              nodes), collapse = ' '), side = 3, line = 1)
        }
      }
    } else {
      if(is.null(nodes) && !is.null(taxon)){
        plot(shape_file, axes = TRUE, main = 'Minimum spanning tree(s) of the terminal node(s)',
             sub = paste0(c('resolution: ', resol[1], ' x ', resol[2]), collapse = ''), cex.main = 0.9,
             cex.sub = 0.7, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
        abline(h = 0, col = 'red', lty = 2) # equator
      } else if(!is.null(nodes) || !is.null(taxon)){
        plot(shape_file, axes = TRUE, main = paste0(c('Mapping a minimum spanning tree of the node(s):',
            nodes), collapse = ' '), sub = paste0(c('resolution: ', resol[1], ' x ', resol[2]), collapse = ''),
            cex.main = 0.9, cex.sub = 0.7, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
        abline(h = 0, col = 'red', lty = 2) # equator
      } else if(is.null(nodes) && is.null(taxon)){
        plot(shape_file, axes = TRUE, main = 'Mapping a minimum spanning tree of the terminal nodes!',
          sub = paste0(c('resolution: ', resol[1], ' x ', resol[2]), collapse = ''),
          cex.main = 0.9, cex.sub = 0.7, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
        abline(h = 0, col = 'red', lty = 2) # equator
      }
    }
  }

  
  # preparing species' msts...
  #######
  ### vetor de acúmulo de espécimes ###
  species <- unique(rownames(coordin))
  qde <- 0
  for(i in species){qde[i] <- length(which(rownames(coordin) == i))}
  qde <- cumsum(qde)
  tabela <- matrix(NA, nrow = dim(coordin), nc = 4)
  #######
  
  if(!is.null(tree)){
    print('5) calculating mst...')
  } else {
    print('4) calculating mst...')
  }
  if(is.null(nodes)){
    
    if(!is.null(taxon)){
      
      if(mintreeall == TRUE){

        for(j in taxon){
          # sppp <- which(unique(rownames(coords)) == taxon[j])
          
          ##### Se quisermos fazer MST de todos os espécimes JUNTOS das espécies escolhidas... 
          idx <- which(names(qde) == j)
          tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- coordin[which(rownames(coordin) == j), 1]
          tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- coordin[which(rownames(coordin) == j), 2]
          tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))
        }
        frame <- as.data.frame(na.exclude(tabela))
        colnames(frame) <- c('species', colnames(coordin[, c(1, 2)]))
        
        Long <- frame[,2]
        Lat <- frame[, 3]
        names(Long) <- frame[,1]
        
        tempo <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(frame), 2, dimnames = list(frame[,1],
                     colnames(frame)[c(2,3)]))
        
        #### shapefile ###
        tempo.d <- as.data.frame(tempo)
        tempo_shape <- lats2Shape(lats = tempo.d)
        # dir.create('out/')
        write.shapefile(tempo_shape, 'out/pointsshape_mintreeall')
        # resul1_shape <- rgdal::readOGR(dsn = 'out/pointsshape_mintreeall.shp', verbose = FALSE)
        # resul1_shape <- readShapeSpatial('tempshape1_out.shp') # shapefile
        resul1_shape <- vect('out/pointsshape_mintreeall.shp', crs = "+proj=longlat +datum=WGS84")
        # proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        # projection(resul1_shape) <- CRS("+proj=longlat +datum=WGS84")
        
        
        ##### MST based on geographic distance #####
        rownames(resul1_shape@coords) <- rownames(tempo.d)
        colnames(resul1_shape@coords) <- c('longitude', 'latitude')
        dista <- earth.dist(lats = resul1_shape)
        mst2 <- dino.mst(dista)
        rownames(mst2) <- rownames(tempo.d)
        colnames(mst2) <- rownames(tempo.d)
        mst_shape <- msn2Shape(msn = mst2, lats = resul1_shape, dist = NULL)
        write.shapefile(mst_shape, 'out/mst_mintreeall')
        
        if(!is.null(tree)){
          print('5) calculating mst... Done')
        } else {
          print('4) calculating mst... Done')
        }
        
        # plotting shapes
        pontos_linha <- shapefile('out/mst_mintreeall.shp', warnPRJ = FALSE)
        proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
        projection(pontos_linha) <- CRS("+proj=longlat +datum=WGS84")
        
        plot(pontos_linha, col = cols1[1], lwd = 3, add = T)
        plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[1:length(taxon)], add = T)
        
        # minimum convex polygon
        if(pol == TRUE){
          conv <- convexhull.xy(tempo)
          plot(conv, add = T, col = adjustcolor(cols1[taxon[1]], transp))
          write.shapefile(conv, 'out/mcp_mintreeall')
        }
        
        # legend
        if(caption == TRUE){
          x <- taxon
          if(sobrepo == TRUE){
            legend(x = extent(shape_file)[2] - 30, y = extent(shape_file)[4], legend = x, pch = 19, col = cols1[x], bty = 'n',
                 pt.cex = 0.6, cex = 0.6, title = 'Adding species...', title.col = 'red')
          } else {
            legend(x = "bottomleft", legend = x, pch = 19, col = cols1[x], title = 'Members of the minimum spanning tree',
             title.col = 'red', pt.cex = 0.6, cex = 0.6)
          }
        }

        # back-transforming lines in points:
        pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular')
        
        # presence-absence matrix:
        ncellras <- ncell(r)
        coor.l <- matrix(NA, nr = ncellras, nc = length(taxon), dimnames = list(seq(1:ncellras),
                      unique(rownames(tempo))[c(1:length(taxon))]))  # tabela com todas as células
        linhasRaster <- rasterize(pontos_linha2, r, field = 1) # raster com as presenças

        writeRaster(linhasRaster, "out/presence_mintreeall.tif",
                    overwrite = TRUE)
    
        plot(linhasRaster, axes = FALSE, legend = FALSE, add = TRUE, col = cols1[1],
             alpha = transp)
        
        # preparing the matrix
        for(i in 1:ncellras){
          if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == FALSE){
            coor.l[i, ] <- 1
          } else if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == TRUE){
            coor.l[i, ] <- 0
          }
        }
        
        coor.l <- na.exclude(coor.l)
        coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
        rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
        write.table(x = coor.ll, file = 'out/pres_abs.txt', sep = '\t')

        return(coor.ll)
        
        
      } else if (mintreeall == FALSE){

        conta <- 0
        coor.l <- matrix(NA, nr = ncellras, nc = length(taxon), dimnames = list(seq(1:ncellras),
                              taxon)[c(1:2)])  # tabela com todas as células
        lista_r <- list()
        
        for(j in taxon){
          conta <- conta + 1
          tempo <- subset(coordin[, 1:2], rownames(coordin) == j)
          
          #### shapefile ###
          tempo.d <- as.data.frame(tempo)
          tempo_shape <- lats2Shape(lats = tempo.d)
          # dir.create('out/')
          write.shapefile(tempo_shape, paste0(c('out/pointshape_', j), collapse = ''))
          resul1_shape <- vect(paste0(c('out/pointshape_', j, '.shp'),
                collapse = ''), crs = "+proj=longlat +datum=WGS84")
          # resul1_shape <- readShapeSpatial('tempshape1_out.shp') # shapefile
          # proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
          
          ##### MST based on geographic distance #####
          rownames(resul1_shape@coords) <- rownames(tempo.d)
          colnames(resul1_shape@coords) <- c('longitude', 'latitude')
          dista <- earth.dist(lats = resul1_shape)
          mst2 <- dino.mst(dista)
          rownames(mst2) <- rownames(tempo.d)
          colnames(mst2) <- rownames(tempo.d)
          mst_shape <- msn2Shape(msn = mst2, lats = resul1_shape, dist = NULL)
          write.shapefile(mst_shape, paste0(c('out/mst_', j), collapse = ''))
          
          if(!is.null(tree)){
            print(paste0(c(conta + 5, ') calculating mst... Done'), collapse = ''))
          } else {
            print(paste0(c(conta + 4, ') calculating mst... Done'), collapse = ''))
          }
          
          # plotting shapes
          pontos_linha <- shapefile(paste0(c('out/mst_', j, '.shp'), collapse = ''),
                                    warnPRJ = FALSE)
          proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
          plot(pontos_linha, col = cols1[j], lwd = 3, add = T)
          plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[j], add = T)
          
          # minimum convex polygon
          if(pol == TRUE){
            conv <- convexhull.xy(tempo)
            write.shapefile(conv, paste0(c('out/mcp_', j), collapse = ''))
            plot(conv, add = T, col = adjustcolor(cols1[j], transp))
          }
          
          idx <- which(names(qde) == j)
          tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- coordin[which(rownames(coordin) == j), 1]
          tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- coordin[which(rownames(coordin) == j), 2]
          tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))
          
          
          # back-transforming lines in points:
          pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular')
          
          
          lista_r[[conta]] <- rasterize(pontos_linha2, r, field = 1) # raster com as presenças
          writeRaster(lista_r[[conta]], paste0(c("out/presence_mst_", j, ".tif"), 
            collapse = ''), format = "GTiff", overwrite = TRUE)
          plot(lista_r[[conta]], axes = FALSE, legend = FALSE, add = TRUE,
               col = cols1[j], alpha = transp)
          # colocando labels:
          teste <- tabelao[j,j] # posições dos táxons do nó
          text(lista_r[[conta]], labels = rep(teste, dim(tempo)[1]), cex = 0.8, pos = 2, col = cols1[j])
          
          # presence-absence matrix:
          # preparing the matrix
          for(i in 1:ncellras){
            if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == FALSE){
              coor.l[i, conta] <- 1
            } else if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == TRUE){
              coor.l[i, conta] <- 0
            }
          }
        }
        
        if(caption == TRUE){
          x <- taxon
          if(sobrepo == TRUE){
            legend(x = extent(shape_file)[2] - 30, y = extent(shape_file)[4], legend = x, pch = 19, col = cols1[x],
                 pt.cex = 0.6, cex = 0.6, title = 'Adding species...', title.col = 'red')
          } else {
              legend(x = "bottomleft", legend = x, pch = 19, col = cols1[x], title = 'Members of the minimum spanning tree(s)',
                title.col = 'red', pt.cex = 0.6, cex = 0.6)
          }
        }
        
        coor.l <- na.exclude(coor.l)
        coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
        rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
        write.table(x = coor.ll, file = 'out/pres_abs.txt', sep = '\t')

        return(coor.ll)
        
      }
    } else if(is.null(taxon)){ # null taxa

      tempo <- coordin
      
      # tempo <- as.data.frame(na.exclude(tempo))
      # colnames(tempo) <- c('species', colnames(coordin[, c(1, 2)]))
      
      Long <- tempo[,1]
      Lat <- tempo[, 2]
      names(Long) <- rownames(tempo)
      
      tempo <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(tempo), 2, dimnames = list(rownames(tempo),
                                  colnames(tempo)[c(1,2)]))
      

      if(mintreeall == TRUE){
        
        #### shapefile ###
        tempo.d <- as.data.frame(tempo)
        tempo_shape <- lats2Shape(lats = tempo.d)
        write.shapefile(tempo_shape, 'out/pointsshape_mintreeall')
        resul1_shape <- vect('out/pointsshape_mintreeall.shp', crs = "+proj=longlat +datum=WGS84")
        # resul1_shape <- readShapeSpatial('tempshape1_out.shp') # shapefile
        proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        
        ##### MST based on geographic distance #####
        rownames(resul1_shape@coords) <- rownames(tempo.d)
        colnames(resul1_shape@coords) <- c('longitude', 'latitude')
        dista <- earth.dist(lats = resul1_shape)
        mst2 <- dino.mst(dista)
        rownames(mst2) <- rownames(tempo.d)
        colnames(mst2) <- rownames(tempo.d)
        mst_shape <- msn2Shape(msn = mst2, lats = resul1_shape, dist = NULL)
        write.shapefile(mst_shape, 'out/mst_mintreeall')
        
        if(!is.null(tree)){
          print('5) calculating mst... Done')
        } else {
          print('4) calculating mst... Done')
        }
        
        # plotting shapes
        pontos_linha <- shapefile('out/mst_mintreeall.shp', warnPRJ = FALSE)
        proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
        plot(pontos_linha, col = cols1[unique(rownames(tempo))[1]], lwd = 3, add = T)
        plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[unique(rownames(tempo))[1]], add = T)
        
        # minimum convex polygon
        if(pol == TRUE){
          conv <- convexhull.xy(tempo)
          plot(conv, add = T, col = adjustcolor(cols1[unique(rownames(tempo))[1]], transp))
          write.shapefile(conv, 'out/mcp_mintreeall')
        }
        
        # legend
        if(caption == TRUE){
          x <- 'All species'
          legend(x = "bottomleft", legend = x, pch = 19, col = cols1[unique(rownames(tempo))[1]],
                 bty = 'o', pt.cex = 0.6, cex = 0.6)
        }
        
        # back-transforming lines in points:
        pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular')
        
        # presence-absence matrix:
        # ncellras <- ncell(r)
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo))), dimnames = list(seq(1:ncellras),
                      unique(rownames(tempo)))[c(1:2)])
        linhasRaster <- rasterize(pontos_linha2, r, field = 1) # raster com as presenças
        writeRaster(linhasRaster, "out/presence_mintreeall.tif", format = "GTiff",
                    overwrite = TRUE)
        plot(linhasRaster, axes = FALSE, legend = FALSE, add = TRUE, col = cols1[1], alpha = transp)
        
        # preparing the matrix
        for(i in 1:ncellras){
          if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == FALSE){
            coor.l[i, ] <- 1
          } else if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == TRUE){
            coor.l[i, ] <- 0
          }
        }
        
        coor.l <- na.exclude(coor.l)
        coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
        rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
        write.table(x = coor.ll, file = 'out/pres_abs.txt', sep = '\t')

        return(coor.ll)
        

      } else if (mintreeall == FALSE){

        conta <- 0
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo))),
                         dimnames = list(seq(1:ncellras), unique(rownames(tempo)))
                                         [c(1:2)])
        lista_r <- list()
        
        for(j in unique(rownames(tempo))){
          conta <- conta + 1
          tempoo <- subset(tempo[, 1:2], rownames(tempo) == j)
          
          #### shapefile ###
          tempo.d <- as.data.frame(tempoo)
          tempo_shape <- lats2Shape(lats = tempo.d)
          write.shapefile(tempo_shape, paste0(c('out/pointshape_', j), collapse = ''))
          resul1_shape <- vect(paste0(c('out/pointshape_', j, '.shp'),
                      collapse = ''), crs = "+proj=longlat +datum=WGS84")
          # resul1_shape <- readShapeSpatial('tempshape1_out.shp') # shapefile
          # proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
          
          # resul1_shape <- vect(resul1_shape, crs = "+proj=longlat")
          
          ##### MST based on geographic distance #####
          rownames(resul1_shape@coords) <- rownames(tempo.d)
          colnames(resul1_shape@coords) <- c('longitude', 'latitude')
          dista <- earth.dist(lats = resul1_shape)
          mst2 <- dino.mst(dista)
          rownames(mst2) <- rownames(tempo.d)
          colnames(mst2) <- rownames(tempo.d)
          mst_shape <- msn2Shape(msn = mst2, lats = resul1_shape, dist = NULL)
          write.shapefile(mst_shape, paste0(c('out/mst_', j), collapse = ''))
          
          if(!is.null(tree)){
            print(paste0(c(conta + 5, ') calculating mst... Done'), collapse = ''))
          } else {
            print(paste0(c(conta + 4, ') calculating mst... Done'), collapse = ''))
          }
          
          # plotting shapes
          pontos_linha <- shapefile(paste0(c('out/mst_', j, '.shp'), collapse = ''),
                                    warnPRJ = FALSE)
          proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
          plot(pontos_linha, col = cols1[j], lwd = 3, add = T)
          plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[j], add = T)
          
          # minimum convex polygon
          if(pol == TRUE){
            conv <- convexhull.xy(tempoo)
            write.shapefile(conv, paste0(c('out/mcp_', j), collapse = ''))
            plot(conv, add = T, col = adjustcolor(cols1[j], transp))
          }
          
          idx <- which(names(qde) == j)
          tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- tempo[which(rownames(tempo) == j), 1]
          tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- tempo[which(rownames(tempo) == j), 2]
          tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))
          
          
          # back-transforming lines in points:
          pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular')
          
          
          lista_r[[conta]] <- rasterize(pontos_linha2, r, field = 1) # raster com as presenças
          plot(lista_r[[conta]], axes = FALSE, legend = FALSE, add = TRUE, col = cols1[j],
               alpha = transp)
          writeRaster(lista_r[[conta]], paste0(c("out/presence_mst_", j, ".tif"), 
            collapse = ''), format = "GTiff", overwrite = TRUE)
  
          # teste <- tabelao[j,j] # posições dos táxons do nó
          teste <- conta
          text(tempoo, labels = rep(teste, dim(tempo)[1]), cex = 0.5, pos = 2, col = cols1[j])
          
          # presence-absence matrix:
          # preparing the matrix
          for(i in 1:ncellras){
            if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == FALSE){
              coor.l[i, conta] <- 1
            } else if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == TRUE){
              coor.l[i, conta] <- 0
            }
          }
        }
        
        # legend
        if(caption == TRUE){
          x <- unique(rownames(tempo))
          legend(x = "bottomleft", legend = x, pch = 19, col = cols1[unique(rownames(tempo))
              [1:length(unique(rownames(tempo)))]], bty = 'o', pt.cex = 0.6, cex = 0.6)
        }
        
        coor.l <- na.exclude(coor.l)
        coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
        rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
        write.table(x = coor.ll, file = 'out/pres_abs.txt', sep = '\t')

        return(coor.ll)
        

      }
    }
    ########### Choosing a node or nodes... ############
  } else if(!is.null(nodes)){
    
    anc <- nodes
    no_num <- 0

    ## looping the nodes:
    for(k in anc){
      no_num <- no_num + 1

      lis <- arv$tip.label[arv$edge[which(arv$edge[,1] == k), 2]]
        
      # criando a condição de checagem dos táxons terminais:
      if (any(is.na(lis)) == FALSE){
        for(j in lis){
          # sppp <- which(unique(rownames(coords)) == taxon[j])
          
          ##### Se quisermos fazer MST de todos os espécimes JUNTOS das espécies escolhidas... 
          idx <- which(names(qde) == j)
          tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- coordin[which(rownames(coordin) == j), 1]
          tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- coordin[which(rownames(coordin) == j), 2]
          tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))
          tabela[(qde[idx - 1] + 1) : qde[idx], 4] <- rep(no_num, (qde[idx] - qde[idx - 1]))
        }

      } else if(any(is.na(lis)) == TRUE){

        # tabelao <- mrca(phy = tree, full = F)
        teste <- which(tabelao %in% anc[no_num]) # posições dos táxons do nó
        tabelao.v <- as.vector(tabelao)
        names(tabelao.v) <- rep(colnames(tabelao), dim(tabelao)[1])
        clado <- unique(names(tabelao.v)[teste])
        
        for(j in clado){
          # sppp <- which(unique(rownames(coords)) == taxon[j])
          
          ##### Se quisermos fazer MST de todos os espécimes JUNTOS das espécies escolhidas... 
          idx <- which(names(qde) == j)
          tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- coordin[which(rownames(coordin) == j), 1]
          tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- coordin[which(rownames(coordin) == j), 2]
          tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))
          tabela[(qde[idx - 1] + 1) : qde[idx], 4] <- rep(no_num, (qde[idx] - qde[idx - 1]))
        }
      }
    }

    # creating a table just for those terminal species  
    if(!is.null(taxon)){
      # preparing species' msts...
      #######
      ### vetor de acúmulo de espécimes ###
      qde_sp <- 0
      for(b in species){qde_sp[b] <- length(which(rownames(coordin) == b))}
      qde_sp <- cumsum(qde_sp)
      tab_sp <- matrix(NA, nrow = dim(coordin), nc = 3)
      #######

      for(z in taxon){
        # sppp <- which(unique(rownames(coordin)) == taxon[j])
        
        ##### Se quisermos fazer MST de todos os espécimes JUNTOS das espécies escolhidas... 
        idx <- which(names(qde_sp) == z)
        tab_sp[(qde_sp[idx - 1] + 1) : qde_sp[idx], 2] <- coordin[which(rownames(coordin) == z), 1]
        tab_sp[(qde_sp[idx - 1] + 1) : qde_sp[idx], 3] <- coordin[which(rownames(coordin) == z), 2]
        tab_sp[(qde_sp[idx - 1] + 1) : qde_sp[idx], 1] <- rep(z, (qde_sp[idx] - qde_sp[idx - 1]))
      }
      tabela <- tabela
    } else if(is.null(taxon)){
      tabela <- tabela
    }

    if(mintreeall == TRUE){
      
      # internal nodes:
      print('adding internal nodes')
      frame.n <- as.data.frame(na.exclude(tabela))
      colnames(frame.n) <- c('species', colnames(coordin[, c(1, 2)]))
      
      Long_n <- frame.n[,2]
      Lat_n <- frame.n[, 3]
      names(Long_n) <- frame.n[,1]
      
      tempo.temp2 <- matrix(cbind(as.numeric(Long_n), as.numeric(Lat_n)), nrow(frame.n), 2,
                           dimnames = list(frame.n[,1], colnames(frame.n)[c(2,3)]))
      #######
      
      if(!is.null(taxon)){

        # terminal nodes:
        print('adding terminal nodes...')
        frame <- as.data.frame(na.exclude(tab_sp))
        colnames(frame) <- c('species', colnames(coordin[, c(1, 2)]))
        
        Long <- frame[,2]
        Lat <- frame[, 3]
        names(Long) <- frame[,1]
        
        tempo.temp <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(frame), 2,
                        dimnames = list(frame[,1], colnames(frame)[c(2,3)]))
        #######
        
        # putting in there alltogether:
        total <- rbind(frame.n, frame)
      
        Long <- total[,2]
        Lat <- total[, 3]
        names(Long) <- total[,1]
        
        tempo <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(total), 2,
                        dimnames = list(total[,1], colnames(total)[c(2,3)]))
      
        #### shapefiles ###
        ## terminal + nodes
        tempo.d <- as.data.frame(tempo)
        tempo_shape <- lats2Shape(lats = tempo.d)
        # dir.create('temp/')
        write.shapefile(tempo_shape, 'out/ancterminal_points_mintreeall')
        resul1_shape <- vect('out/ancterminal_points_mintreeall.shp',
         crs = "+proj=longlat +datum=WGS84")
        # proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        # projection(resul1_shape) <- CRS("+proj=longlat +datum=WGS84")
        
        #only terminals:
        tempo_temp.d <- as.data.frame(tempo.temp)
        tempo_shape.temp <- lats2Shape(lats = tempo_temp.d)
        write.shapefile(tempo_shape.temp, 'out/points_mintreeall_onlyterminal')
        resul1_shape.temp <- vect('out/points_mintreeall_onlyterminal.shp',
          crs = "+proj=longlat +datum=WGS84")
        # proj4string(resul1_shape.temp) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        # projection(resul1_shape.temp) <- CRS("+proj=longlat +datum=WGS84")
      }
      
      #only internal nodes:
      tempo_temp.d2 <- as.data.frame(tempo.temp2)
      tempo_shape.temp <- lats2Shape(lats = tempo_temp.d2)
      write.shapefile(tempo_shape.temp, 'out/points_mintreeall_onlyinternal')
      resul2_shape.temp <- vect('out/points_mintreeall_onlyinternal.shp',
       crs = "+proj=longlat +datum=WGS84")
      # proj4string(resul2_shape.temp) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
      # projection(resul2_shape.temp) <- CRS("+proj=longlat +datum=WGS84")
      
      if(!is.null(taxon)){
        ##### total MST based on geographic distance #####
        rownames(resul1_shape@coords) <- rownames(tempo.d)
        colnames(resul1_shape@coords) <- c('longitude', 'latitude')
        dista <- earth.dist(lats = resul1_shape)
        mst2 <- dino.mst(dista)
        rownames(mst2) <- rownames(tempo.d)
        colnames(mst2) <- rownames(tempo.d)
        mst_shape <- msn2Shape(msn = mst2, lats = resul1_shape, dist = NULL)
        write.shapefile(mst_shape, 'out/mst_ancterminal_mintreeall')
        
        if(!is.null(tree)){
          print('5) calculating total mst... Done')
        } else {
          print('4) calculating total mst... Done')
        }
        
        # plotting shapes
        pontos_linha <- shapefile('out/mst_ancterminal_mintreeall.shp', warnPRJ = FALSE)
        proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
        projection(pontos_linha) <- CRS("+proj=longlat +datum=WGS84")
        
        plot(pontos_linha, col = 'red', lwd = 3, lty = 2, add = T)
        plot(resul2_shape.temp, cex = 1.5, pch = 'o', add = T)
        # colocando labels:
        text(resul2_shape.temp, labels = rep(nodes, dim(tempo)[1]), cex = 0.8, pos = 1,
             col = cols1[lis])
        
        # plotting only terminal:
        plot(resul1_shape.temp, cex = 1.5, pch = 21, bg = cols1[taxon], add = T)
        for(nomes in taxon){
          teste <- tabelao[nomes,nomes] # posições dos táxons do nó
          text(tempo_temp.d, labels = rep(teste, dim(tempo_temp.d)[1]), cex = 0.8, pos = 2, col = cols1[nomes])
        }
        
        
        # minimum convex polygon
        if(pol == TRUE){
          conv <- convexhull.xy(tempo)
          plot(conv, add = T, col = adjustcolor('gray', transp))
          write.shapefile(conv, 'out/mcp_mintreeall_ancterminal')
        }
        
        # legend
        if(caption == TRUE){
          x <- unique(rownames(tempo))
          xx <- unique(rownames(tempo.temp))
          if(length(anc) == 1){
            if(sobrepo == TRUE){
              legend(x = extent(shape_file)[2] - 30, y = extent(shape_file)[4], legend = x,
                 pch = 'o', col = cols1[lis], title = paste0(c('Adding descendents from node ',
                   nodes, '...'), collapse = ''), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = extent(shape_file)[1], y = extent(shape_file)[4], legend = xx,
                     pch = 'o', col = cols1[xx], title = 'Adding terminal(s)...', 
                     title.col = 'blue', pt.cex = 0.6, cex = 0.6)
            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x,
                 pch = 'o', col = cols1[lis], title = paste0(c('Descendents from node ', nodes),
                        collapse = ''), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = extent(shape_file)[1], y = extent(shape_file)[4], legend = xx,
                     pch = 'o', col = cols1[xx], title = 'Adding terminal(s)...', 
                     title.col = 'blue', pt.cex = 0.6, cex = 0.6)
            }
          } else if(length(anc) > 1){
            if(sobrepo == TRUE){
              legend(x = extent(shape_file)[2] - 30, y = extent(shape_file)[4], legend = x,
                     pch = 'o', col = cols1[lis], title = paste0(c('Adding descendents from nodes',
                nodes, '...'), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = extent(shape_file)[1], y = extent(shape_file)[4], legend = xx,
                     pch = 'o', col = cols1[xx], title = 'Adding terminal(s)...', 
                     title.col = 'blue', pt.cex = 0.6, cex = 0.6)
            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x,
                    pch = 'o', col = cols1[lis], title = paste0(c('Descendents from nodes', nodes),
                    collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = extent(shape_file)[1], y = extent(shape_file)[4], legend = xx,
                     pch = 'o', col = cols1[xx], title = 'Adding terminal(s)...', 
                     title.col = 'blue', pt.cex = 0.6, cex = 0.6)
            }
          }
        }

        # back-transforming lines in points:
        pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular')
        
        # presence-absence matrix:
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo))), dimnames = list(seq(1:ncellras),
                  unique(rownames(tempo)))[c(1:2)])  # tabela com todas as células
        linhasRaster <- rasterize(pontos_linha2, r, field = 1) # raster com as presenças
        writeRaster(linhasRaster, "out/presence_mst_mintreeall_ancterminal.tif",
                    format = "GTiff", overwrite = TRUE)
        # plot(linhasRaster, axes = FALSE, legend = FALSE, add = TRUE, col = cols1[1],
        #     alpha = transp)
        
        # preparing the matrix
        for(i in 1:ncellras){
          if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == FALSE){
            coor.l[i, ] <- 1
          } else if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == TRUE){
            coor.l[i, ] <- 0
          }
        }
      } else {
        ##### MST of the internal nodes based on geographic distance #####
        rownames(resul2_shape.temp@coords) <- rownames(tempo_temp.d2)
        colnames(resul2_shape.temp@coords) <- c('longitude', 'latitude')
        dista <- earth.dist(lats = resul2_shape.temp)
        mst2 <- dino.mst(dista)
        rownames(mst2) <- rownames(tempo_temp.d2)
        colnames(mst2) <- rownames(tempo_temp.d2)
        mst_shape <- msn2Shape(msn = mst2, lats = resul2_shape.temp, dist = NULL)
        write.shapefile(mst_shape, 'out/mst_mintreeall_onlyinternal')
        
        if(!is.null(tree)){
          print('5) calculating mst from internal node(s)... Done')
        } else {
          print('4) calculating mst from internal node(s)... Done')
        }
        
        # plotting shapes
        pontos_linha <- shapefile('out/mst_mintreeall_onlyinternal.shp', warnPRJ = FALSE)
        proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
        projection(pontos_linha) <- CRS("+proj=longlat +datum=WGS84")
        
        plot(pontos_linha, col = 'red', lwd = 2, lty = 2, add = T)
        plot(resul2_shape.temp, cex = 1.5, pch = 'o', add = T)
        # colocando labels:
        text(resul2_shape.temp, labels = rep(nodes, dim(tempo_temp.d2)[1]), cex = 0.8, pos = 1,
             col = cols1[lis])
                
        
        # minimum convex polygon
        if(pol == TRUE){
          conv <- convexhull.xy(tempo_temp.d2)
          plot(conv, add = T, col = adjustcolor('gray', transp))
          write.shapefile(conv, 'out/mcp_mintreeall_ancterminal')
        }


        # legend
        if(caption == TRUE){
          x <- unique(rownames(tempo.temp2))
          if(length(anc) == 1){
            if(sobrepo == TRUE){
              legend(x = extent(shape_file)[2] - 30, y = extent(shape_file)[4], legend = x,
                 pch = 'o', title = paste0(c('Adding descendents from node',
                   nodes, '...'), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)

            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x,
                pch = 'o', title = paste0(c('Descendents from node', nodes),
                collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
      
            }
          } else if(length(anc) > 1){
            if(sobrepo == TRUE){
              legend(x = extent(shape_file)[2] - 30, y = extent(shape_file)[4], legend = x,
                    pch = 'o', title = paste0(c('Adding descendents from nodes',
                    nodes, '...'), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)

            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x,
                  pch = 'o', title = paste0(c('Descendents from nodes', nodes),
                  collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
      
            }
          }
        }

        # back-transforming lines in points:
        pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular')
        
        # presence-absence matrix:
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo.temp2))), dimnames = list(seq(1:ncellras),
                  unique(rownames(tempo.temp2)))[c(1:2)])  # tabela com todas as células
        linhasRaster <- rasterize(pontos_linha2, r, field = 1) # raster com as presenças
        writeRaster(linhasRaster, "out/presence_mst_mintreeall_ancterminal.tif",
                    format = "GTiff", overwrite = TRUE)
        # plot(linhasRaster, axes = FALSE, legend = FALSE, add = TRUE, col = cols1[1],
        #     alpha = transp)
        
        # preparing the matrix
        for(i in 1:ncellras){
          if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == FALSE){
            coor.l[i, ] <- 1
          } else if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == TRUE){
            coor.l[i, ] <- 0
          }
        }
      }

      coor.l <- na.exclude(coor.l)
      coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
      rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
      write.table(x = coor.ll, file = 'out/pres_abs.txt', sep = '\t')

      return(coor.ll)
      

    } else if(mintreeall == FALSE){
      # internal nodes:
      print('adding internal nodes...')
      frame.n <- as.data.frame(na.exclude(tabela))
      # frame.n <- as.data.frame(tabela)
      colnames(frame.n) <- c('species', colnames(coordin[, c(1, 2)]), 'nodes')
      
      Long_n <- frame.n[,2]
      Lat_n <- frame.n[, 3]
      nodd <- frame.n[, 4]
      names(Long_n) <- frame.n[,1]
      #######

      # terminal nodes:
      if(!is.null(taxon)){
        print('adding terminal nodes...')
        frame <- as.data.frame(na.exclude(tab_sp))
        colnames(frame) <- c('species', colnames(coordin[, c(1, 2)]))
        
        Long <- frame[,2]
        Lat <- frame[, 3]
        names(Long) <- frame[,1]
      }
        #######

      # internal nodes:
      tempo.n <- matrix(cbind(as.numeric(Long_n), as.numeric(Lat_n), as.numeric(nodd)), nrow(frame.n), 3,
                   dimnames = list(frame.n[,1], colnames(frame.n)[c(2:4)]))

      if(!is.null(taxon)){
        #terminal nodes:
        tempo <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(frame), 2,
                        dimnames = list(frame[,1], colnames(frame)[c(2,3)]))
      }

      if(!is.null(taxon)){
        #let's put in there alltogether!
        coor.l <- matrix(NA, nr = ncellras, nc = sum(length(unique(rownames(tempo.n))), length(unique(rownames(tempo)))),
         dimnames = list(seq(1:ncellras), c(unique(rownames(tempo.n)), unique(rownames(tempo))))[c(1:2)])  # tabela com todas as células
      } else {
        #let's put in there alltogether!
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo.n))),
         dimnames = list(seq(1:ncellras), unique(rownames(tempo.n)))[c(1:2)])  # tabela com todas as células

      }
      
      # first, to internal nodes...
      # lista_r <- list()
      conta <- 0
      
      tempoo <- tempo.n[, 1:2]
      
      #### shapefile ###
      tempo.d <- as.data.frame(tempoo)
      tempo_shape <- lats2Shape(lats = tempo.d)
      # dir.create('out/')
      write.shapefile(tempo_shape, paste0(c('out/pointshape_node', nodes), collapse = ''))
      resul1_shape <- vect(paste(c('out/pointshape_node', nodes, '.shp'),
                                 collapse = ''), crs = "+proj=longlat +datum=WGS84")
      # resul1_shape <- readShapeSpatial('tempshape1_out.shp') # shapefile
      # proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
      # suppressWarnings(projection(resul1_shape) <- CRS("+proj=longlat +datum=WGS84"))
      
      ##### MST based on geographic distance #####
      # rownames(resul1_shape@coords) <- rownames(tempo.d)
      # colnames(resul1_shape@coords) <- c('longitude', 'latitude')
      dista <- earth.dist(lats = as(resul1_shape, 'Spatial'))
      mst2 <- dino.mst(dista)
      rownames(mst2) <- rownames(tempo.d)
      colnames(mst2) <- rownames(tempo.d)
      lats <- cbind(resul1_shape$LONG,
                    resul1_shape$LAT)
      rownames(lats) <- rownames(tempo.d)
      colnames(lats) <- c('longitude', 'latitude')
      mst_shape <- msn2Shape(msn = mst2, lats = lats, dist = NULL)
      write.shapefile(mst_shape, paste0(c('out/mst_node_', nodes), collapse = ''))
      
      if(!is.null(tree)){
        print(paste0(c(conta + 5, ') calculating mst... Done'), collapse = ''))
        contass <- conta + 5
      } else {
        print(paste0(c(conta + 4, ') calculating mst... Done'), collapse = ''))
        contass <- conta + 4
      }
      
      # plotting shapes
      pontos_linha <- vect(paste0(c('out/mst_node_', nodes, '.shp'), collapse = ''),
                           crs = "+proj=longlat +datum=WGS84")
      # proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
      plot(pontos_linha, col = 'red', lwd = 3, lty = 2, add = T)
      plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[5], add = T)
      
      # labels:
      for(j in 1:length(nodes)){
        text(tempo.n[which(tempo.n [, 3] == j), c(1, 2)],
             labels = rep(nodes[j], dim(tempo.n[which(tempo.n [, 3] == j), c(1, 2)])[1]),
             cex = 0.6, pos = 1, col = 'black')
      }
      
      # minimum convex polygon
      if(pol == TRUE){
        conv <- convexhull.xy(tempoo)
        plot(conv, add = T, col = adjustcolor(cols1[j], transp))
        write.shapefile(conv, paste0(c('out/mcp_node_', nodes), collapse = ''))
      }              
      
      # back-transforming lines in points:
      suppressWarnings(pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular'))
      
      lista_r <- rasterize(pontos_linha2, r, field = 1) # raster com as presenças
      writeRaster(lista_r, paste0(c("out/presence_mst_", nodes, ".tif"), 
                                  collapse = ''), overwrite = TRUE)
      
      plot(lista_r, axes = FALSE, legend = FALSE, add = TRUE,
           col = cols1[5], alpha = transp)
      
      # presence-absence matrix:
      # preparing the matrix
      for(k in unique(rownames(tempo.n))){
        for(i in 1:ncellras){
          if(is.na(r[i]) == FALSE && is.na(lista_r[i]) == FALSE){
            coor.l[i, k] <- 1
          } else if(is.na(r[i]) == FALSE && is.na(lista_r[i]) == TRUE){
            coor.l[i, k] <- 0
          }
        }
      }
  

      if(!is.null(taxon)){
        #finally, the terminal nodes:   
        contas <- 0
        for(j in taxon){
          contas <- contas + 1
          conta <- conta + 1
          tempoo <- subset(tempo[, 1:2], rownames(tempo) == j)
          
          #### shapefile ###
          tempo.d <- as.data.frame(tempoo)
          tempo_shape <- lats2Shape(lats = tempo.d)
          # dir.create('temp/')
          write.shapefile(tempo_shape, paste0(c('out/pointshape_', j), collapse = ''))
          resul1_shape <- rgdal::readOGR(dsn = paste0(c('out/pointshape_', j, '.shp'),
                                            collapse = ''), verbose = FALSE)
          # resul1_shape <- readShapeSpatial('tempshape1_out.shp') # shapefile
          proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
          projection(resul1_shape) <- CRS("+proj=longlat +datum=WGS84")
          
          ##### MST based on geographic distance #####
          rownames(resul1_shape@coords) <- rownames(tempo.d)
          colnames(resul1_shape@coords) <- c('longitude', 'latitude')
          dista <- earth.dist(lats = resul1_shape)
          mst2 <- dino.mst(dista)
          rownames(mst2) <- rownames(tempo.d)
          colnames(mst2) <- rownames(tempo.d)
          mst_shape <- msn2Shape(msn = mst2, lats = resul1_shape, dist = NULL)
          write.shapefile(mst_shape, paste0(c('out/mst_', j), collapse = ''))
          
          print(paste0(c(contas + contass, ') calculating mst... Done'), collapse = ''))
          
          # plotting shapes
          pontos_linha <- shapefile(paste0(c('out/mst_', j, '.shp'), collapse = ''),
                                    warnPRJ = FALSE)
          proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
          plot(pontos_linha, col = cols1[j], lwd = 3, lty = 2, add = T)
          # plotting only terminal:
          plot(resul1_shape, cex = 1.5, pch = 21, bg = cols1[j], add = T)
          
          teste <- tabelao[j, j] # posições dos táxons do nó
          text(tempoo, labels = rep(teste, dim(tempoo)[1]), cex = 0.8, pos = 2, col = cols1[j])

          # minimum convex polygon
          if(pol == TRUE){
            conv <- convexhull.xy(tempoo)
            plot(conv, add = T, col = adjustcolor(cols1[j], transp))
            write.shapefile(conv, paste0(c('out/mcp_', j), collapse = ''))
          }              
                
          # back-transforming lines in points:
          pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular')
                  
          lista_r[[conta]] <- rasterize(pontos_linha2, r, field = 1) # raster com as presenças
          writeRaster(lista_r[[conta]], paste0(c("out/presence_mst_", j, ".tif"), 
            collapse = ''), format = "GTiff", overwrite = TRUE)
          
          plot(lista_r[[conta]], axes = FALSE, legend = FALSE, add = TRUE, col = cols1[j],
               alpha = transp)
                
          # presence-absence matrix:
          # preparing the matrix
          for(i in 1:ncellras){
            if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == FALSE){
              coor.l[i, conta] <- 1
            } else if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == TRUE){
              coor.l[i, conta] <- 0
            }
          }
        }
      }

      if(!is.null(taxon)){
        # legend
        if(caption == TRUE){
          x <- unique(rownames(tempo.n))
          xx <- unique(rownames(tempo))
          xxx <- unique(tempo.n[,3])
          if(length(anc) == 1){
            if(sobrepo == TRUE){
              legend(x = extent(shape_file)[2] - 30, y = extent(shape_file)[4], legend = x,
                     pch = 19, col = cols1[x], title = paste0(c('Adding descendents from the node',
                   nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'bottomright', legend = xx,
                     pch = 19, col = cols1[xx], title = 'Adding terminal node(s)...', title.col = 'blue',
                 pt.cex = 1.2, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)
              
            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x, pch = 19, col = cols1[x], title = paste0(c('Descendents from the node',
               nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'bottomright', legend = xx,
                 pch = 19, col = cols1[xx], title = 'Adding terminal node(s)...',
                   title.col = 'blue', pt.cex = 1.2, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)
              
            }
          } else if(length(anc) > 1){
            if(sobrepo == TRUE){
              legend(x = extent(shape_file)[2] - 30, y = extent(shape_file)[4], legend = x,
                     pch = 19, col = cols1[x], title = paste0(c('Adding descendents from the nodes',
                   nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'bottomright', legend = xx,
                 pch = 19, col = cols1[xx], title = 'Adding terminal node(s)...',
                 title.col = 'blue', pt.cex = 1.2, cex = 0.6, bty = 'n')
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)
            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x, pch = 19, col = cols1[x], title = paste0(c('Descendents from the nodes',
               nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'bottomright', legend = xx,
                 pch = 19, col = cols1[xx], title = 'Adding terminal node(s)...',
                 title.col = 'blue', pt.cex = 1.2, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)
            }
          }
        }
      } else {
        # legend
        if(caption == TRUE){
          x <- unique(rownames(tempo.n))
          xxx <- unique(tempo.n[,3])
          if(length(anc) == 1){
            if(sobrepo == TRUE){
              legend(x = extent(shape_file)[2] - 30, y = extent(shape_file)[4], legend = x,
                  pch = 19, col = cols1[x], title = paste0(c('Adding descendents from the node',
                  nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)
      
            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x, pch = 19, col = cols1[x], title = paste0(c('Descendents from the node',
               nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)

            }
          } else if(length(anc) > 1){
            if(sobrepo == TRUE){
              legend(x = extent(shape_file)[2] - 30, y = extent(shape_file)[4], legend = x,
                 pch = 19, col = cols1[x], title = paste0(c('Adding descendents from the nodes',
                   nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)
      
            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x, pch = 19, col = cols1[x], title = paste0(c('Descendents from the nodes',
               nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)
          
            }
          }
        }

      }       

      coor.l <- na.exclude(coor.l)
      coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
      rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
      write.table(x = coor.ll, file = 'out/pres_abs.txt', sep = '\t')

      return(coor.ll)
      
    }
  }
}