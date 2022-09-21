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
  # plot(cropped_map, xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = T)
  r <- rasterize(shapeFile, mask.raster, fun = 'first')
  proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  r <- merge(r, mask.raster)
  ncellR <- ncell(r)
  r[r > 0] <- NA
  for(i in 1:length(gridCell)){
    for(j in 1:ncellR){
      if(j == gridCell[i]){
        values(r)[j] <- 1
      }
    }
  }
  
  map.r <- as.data.frame(rasterToPoints(r))
  pontosRaster <- rasterize(cbind(map.r$x, map.r$y), r, mask = T)
  # plot(pontosRaster)
  gridNumber <- which(pontosRaster@data@values == 1)
  # map.r$gridNumber <- which(pontosRaster@data@values == 1)
  
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
    points(tempT[, c(2, 3)], pch = 16, col = adjustcolor(col = names(cols1[contando]), alpha = 0.6), cex = 1)
  }
  legend(x = "bottomleft", legend = speciesN, pch = 16, col = names(cols1[1:length(speciesN)]),
         title = 'Species',  title.col = 'red', pt.cex = 0.6, cex = 0.6)
  
  plot(pontosRaster, add = T, axes = F, legend = F, col = colores, alpha = 0.40)
  
  map.r$gridNumber <- which(pontosRaster@data@values == 1)
  
  text(subset(map.r[,c(1, 2)], map.r$layer == 1), labels = gridNumber, cex = (2 * resolut[1]) / resolut[1],
       col = adjustcolor(col = 'black', alpha = 0.6), font = 2)
}