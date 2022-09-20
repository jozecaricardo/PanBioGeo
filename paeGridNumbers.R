gridNumbers <- function(shape_file, resol, coordina, gridCell, transp = 0.8, xmin = NULL,
                        xmax = NULL, ymin = NULL, ymax = NULL){
  
  ######################
  ##### shape file #####
  # shapeFile <- asul
  resolut <- resol
  
  coordin <- matrix(as.matrix(coordina[, c(2, 3)]), nrow(coordina),
                    2, dimnames = list(coordina[,1], colnames(coordina)[c(2, 3)]))
  
  cols1 <- setNames(n = viridis(length(unique(rownames(coordin)))),
                    unique(rownames(coordin)))
  
  xmin = xmin; xmax = xmax; ymin = ymin; ymax = ymax
  
  
  grid <- raster(extent(shapeFile), resolution = resolut, crs = CRS("+proj=longlat +datum=WGS84"))
  grid <- raster::extend(grid, c(1, 1))
  gridPolygon <- rasterToPolygons(grid)
  # suppressWarnings(proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84")) # datum WGS84
  #proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  projection(gridPolygon) <- CRS("+proj=longlat +datum=WGS84")
  
  
  # clipping the intersected cells:
  cropped_map <- raster::intersect(gridPolygon, shapeFile)
  # plot(cropped_map)
  # cropped_map <- vect(cropped_map, geom = c("x", "y"), crs = crs("+proj=longlat +datum=WGS84"))
  
  #########
  # producing a raster of the shapefile
  mask.raster <- raster(extent(gridPolygon), resolution = resolut,
                        crs = CRS("+proj=longlat +datum=WGS84"))
  r <- rasterize(shapeFile, mask.raster)
  proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  # mask.raster[is.na(mask.raster)] <- 0
  r <- merge(r, mask.raster)
  # plot(r)
  
  ncellR <- ncell(r)
  # 
  # # set the cells associated with the shapefile to the specified value
  r[r > 0] <- NA
  
  # data frame com os grids escolhidos
  grid_n <- gridCell
  frameTemp <- data.frame(spp = 'numEscolhido', grid_n = grid_n, row.names = 1:length(grid_n))
  # frameTemp
  
  # matrix with species names and grid numbers:
  resulPaeRaster <- matrix(as.matrix(frameTemp[,2]), nrow(frameTemp),
                           1, dimnames = list(frameTemp[,1], colnames(frameTemp)[2]))
  
  for(j in 1:dim(resulPaeRaster)[1]){
    gTrack <- resulPaeRaster[j, 1]
    values(r)[gTrack] <- 1
  }
  # r[is.na(r)]<-0
  
  #convert the raster to points for plotting the number of a grid
  map.r <- as.data.frame(rasterToPoints(r))
  pontosRaster <- rasterize(cbind(map.r$x, map.r$y), r, field = 1) # raster com as presenças
  
  writeRaster(pontosRaster, 'out/resultPAE_PCE.tif', format = "GTiff", overwrite = TRUE)
  
  colores <- hcl.colors(1)
  brks <- seq(0, 1, by=1)
  
  plot(cropped_map, axes = TRUE, main = 'PAE-PCE final result with numbered grids',
       sub = paste0(c('resolution: ', resolut[1], ' x ', resolut[2]), collapse = ''), cex.main = 0.9,
       cex.sub = 0.7, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  speciesN <- unique(rownames(coordin))
  nspeciesN <- length(speciesN)
  contando <- 0
  for(i in speciesN){
    contando <- contando + 1
    tempT <- subset(coordina, coordina$spp == i)
    points(tempT[, c(2, 3)], pch = 16, col = adjustcolor(col = names(cols1[contando]), alpha = transp), cex = 1)
  }
  legend(x = "bottomleft", legend = speciesN, pch = 16, col = names(cols1[1:length(speciesN)]),
         title = 'Species',  title.col = 'red', pt.cex = 0.6, cex = 0.6)
  
  plot(r, add = T, axes = F, legend = F, col = colores, alpha = transp/2, breaks = brks)
  
  map.r$gridNumber <- which(pontosRaster@data@values == 1)
  
  text(map.r[,c(1, 2)], labels = map.r$gridNumber, cex = (2 * resolut[1]) / resolut[1], col = 'black', font = 2)
  
}


