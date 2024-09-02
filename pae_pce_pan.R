#### PAE-PCE analysis ####
# preabsMat = the presence-absence matrix from terminal_node function
# rou = number of iterations in the parsimony with N runs
pae_pce <- function(preabsMat, shapeFile, resolut, N = NULL,
                    gridView = FALSE, labelGrid = FALSE, legendSpecies = TRUE, sobrepo = FALSE,
                    xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL){
  # install.packages('ctv')
  # install.packages("stringr")
  # install.packages('stringi')
  # install.packages('TreeSearch')
  #library(stringr)
  # library(stringi)
  # library(ctv)
  library(phangorn)
  library(raster)
  library(viridis)
  library(TreeSearch) # TBR search
  # install.views('Phylogenetics')
  
  #Convert the character matrix to a phyDat object 
  # in order to infer a tree. By assigning the matrix 
  # type as 'user' we can specify the components of the 
  # matrix (in this case, binary 1s and 0s).
  
  # matTemp <- preabsMat
  
  tempMatrix <- as.phyDat(preabsMat, type = 'USER', levels = c(0, 1))
  
  colu <- ncol(preabsMat)
  
  # creating a restriction:
  ciVec <- rep(1.0, colu)
  riVec <- rep(1.0, colu)
  
  
  # contagem <- 0
  
  lista <- list()
  listaR <- list()
  
  print('Please, remember that the resulting rasters will be kept in your setted directory!')
  
  if(sobrepo == TRUE){
    print('Please do not forgive that grids are in different colors!')
    if(is.null(N)){
      stop('You should to indicate the number of runs in the N argument!')
    } else {
      print('The N argument is equal to 1 as default!')
    }
  }
  
  ######################
  ##### shape file #####
  shapeFile <- as(shapeFile, 'Spatial')
  grid <- raster(extent(shapeFile), resolution = resolut, crs = CRS("+proj=longlat +datum=WGS84"))
  grid <- raster::extend(grid, c(1, 1))
  gridPolygon <- rasterToPolygons(grid)
  crs(gridPolygon) <- "+proj=longlat +datum=WGS84" # datum WGS84
  #proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  # projection(gridPolygon) <- CRS("+proj=longlat +datum=WGS84")
  
  
  # clipping the intersected cells:
  cropped_map <- raster::intersect(gridPolygon, shapeFile)
  
  # producing a raster of the shapefile
  mask.raster <- raster(extent(shapeFile), resolution = resolut,
                        crs = CRS("+proj=longlat +datum=WGS84"))
  # suppressWarnings(r <- rasterize(shapeFile, mask.raster))
  r <- rasterize(shapeFile, mask.raster)
  proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  # mask.raster[is.na(mask.raster)] <- 0
  r <- merge(r, mask.raster)


  if(sobrepo == FALSE){
    colores <- viridis(n = length(r[!is.na(r)]), option = 'A')
  } else if(sobrepo == TRUE){
    colores <- viridis(n = length(r[!is.na(r)]), option = 'H')
  }
  
  if(sobrepo == FALSE){
    plot(shapeFile, axes = TRUE, main = 'Generalized tracks converted to grids',
         sub = paste0(c('resolution: ', resolut[1], ' x ', resolut[2]), collapse = ''), cex.main = 0.9,
         cex.sub = 0.7, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
    # abline(h = 0, col = 'red', lty = 2) # equator
    
    iterat <- 0
  } else if(sobrepo == TRUE){
    iterat <- 0
  }
  
  # set the cells associated with the shapfile to the specified value
  r[r == 1] <- NA
  
  
  ####### Detecting those grids being supported by synapomorphies ########
  
  #---------------------------------- #### PAE-PCE ###
  # matTemp <- list()
  contagem <- 1
  conta <- 0

  while(length(which((ciVec == 1) == TRUE)) > 1 && dim(preabsMat)[2] > 1){
      
    # We can get a random starting tree for a parsimony search using the
    #  random.addition function:
    # ra.tre <- njs(dist.hamming(tempMatrix))
    
    if(dim(as.data.frame(tempMatrix))[1] < 2){
      stop('This iteration had to be stopped because there is no more synapomorphies!')
    }
    
    ### TBR + RATCHET searches
    dm <- dist.hamming(tempMatrix)
    ra.tre <- NJ(dm)
    
    # Use the pratchet function (parsimony ratchet) to find the most
    # parsimonious tree. Specifying k means the algorithm will search   
    # through k possible trees to find the most parsimonious solution.
    treeIt1 <- list()
    # teste1 <- list()
    
    for(rou in 1:N){
      treeIt1[[rou]] <- optim.parsimony(tree = ra.tre, data = tempMatrix,
                                        method = 'fitch', rearrangements = 'SPR')
      # teste1[[rou]] <- multi2di(treeIt1[[rou]])
      # teste1[[rou]] <- TBR(teste1[[rou]])
      treeIt1[[rou]] <- pratchet(tempMatrix, start = treeIt1[[rou]], maxit=1000,
                              minit=100, k=10, trace=1)
      
      treeIt1[[rou]] <- consensus(treeIt1[[rou]])
      
      ra.tre <- treeIt1[[rou]]
    }
    
    if(is.null(treeIt1)){
      stop('This iteration had to be stopped because there is no more topologies!')
    }
    
    
    
    treeT <- consensus(treeIt1) # 'majority rule': change p to 0.5
    
    # treeComp <- optim.parsimony(tree = treeT, data = tempMatrix, method = 'fitch',
    #                               rearrangements = 'SPR')
    # treeComp <- multi2di(treeComp)
    
    # scoreTreeInit <- parsimony(treeComp, tempMatrix, method = 'fitch')
    
    
    # treeIt <- list()
    # treeIt <- TBR(treeComp)
    # treeIt <- consensus(treeIt1, p = 1.0) # 'majority rule': change p to 0.5
    
    # treeIt2[[times]] <- optim.parsimony(tree = treeT, data = tempMatrix, method = 'fitch',
    #                            rearrangements = 'SPR')
    # parsimony(c(treeIt1, treeIt2), tempMatrix, method = 'fitch')
    
    # Root the tree by the designated outgroup 
    # (write the species name as it appears in the tree, 
    # but add an underscore to fill the space).
    strictTree <- root(treeT, outgroup = "ROOT", resolve.root = TRUE)
    # x11()
    # plot(strictTree)
    # assign edge length
    # treeRatchet <- acctran(tree = treeIt2, data = tempMatrix)
    
    # ancestral states:
    # ancResul <- ancestral.pars(tree = strictTree, data = tempMatrix)
    # plotAnc(tree = strictTree, data = ancResul, attr(ancResul, "index")[2])
    
    # CI of each character (i.e., taxon):
    ciVec <- CI(strictTree, data = tempMatrix, sitewise = TRUE)
    # ciVec[which(is.na(ciVec))] <- 0
    
    # RI of each character (i.e., taxon):
    riVec <- RI(strictTree, data = tempMatrix, sitewise = TRUE)
    # riVec[which(is.na(riVec))] <- 0
    
    if(all(is.na(riVec) == TRUE)){
      print('This analysis had to be stopped because there is no more synapomorphies!')
      break # remaining body of loop will not be executed!
    }
    
    if(any(riVec == 1.0) == FALSE){
      print('This analysis had to be stopped because there is no more synapomorphies!')
      break # remaining body of loop will not be executed!
    }
    
    
    #	colu <- colu - 1
    #}
    
    ## first load packages
    require(phytools)
    tree_temp <- strictTree
    if(is.null(tree_temp$edge.length) == TRUE){
      print('...creating the branch lengths of a tree equal to one...')
      tree_temp <- compute.brlen(tree_temp, 1)
      # tree_temp <- acctran(tree_temp, tempMatrix)
    }
    
    #Discrete mapping
    inter1 <- which(ciVec == 1)
    inter2 <- which(riVec == 1)
    # inter3 <- which(riVec > 0)
    # interT1 <- inter1
    interT1 <- intersect(inter1, inter2)

    if(length(interT1) == 0 && contagem == 1){
      plot(cropped_map, add = TRUE)
      stop('This analysis has no indicative of generalized tracks!')
    }
    
    if(length(interT1) < 2){
      print('This analysis had to be stopped because there is no more synapomorphies!')
      break # remaining body of loop will not be executed!
    }
    
    # coisas que tiveram IC e IR iguais a 1...
    lista_T <- list()
    # inter3 <- which(ciVec == 1); inter4 <- which(riVec == 1)
    # interTT <- intersect(inter3, inter4)
    matTemp2 <- preabsMat[, interT1] # temporary matrix with non-homoplastic synapomorphies
    # matTemp3 <- matTemp2[, -interT2]
    matTemp2 <- as.matrix(matTemp2)
    # matTemp3 <- as.matrix(matTemp3)
    
    nomesCOLsim <- colnames(matTemp2) # species which are synapomorphies
    if(length(nomesCOLsim) < 2){
      print('This analysis had to be stopped because there is no more synapomorphies!')
      break # remaining body of loop will not be executed!
    }
    
    resulT <- dim(matTemp2)[2]
    
    # common grids (synapomorphies)
    roundRes1 <- rownames(matTemp2[rowSums(matTemp2) <= resulT, , drop = F])
    roundRes2 <- rownames(matTemp2[rowSums(matTemp2) > 1, , drop = F])
    roundRes <- intersect(roundRes1, roundRes2)
    if(is.null(roundRes) && contagem == 1){
      plot(cropped_map, add = TRUE)
      stop('This analysis has no indicative of generalized tracks!')
    } else if(is.null(roundRes) && contagem > 1){
      print(paste0(c('The iteration number ', contagem, ' has finished without results.'), collapse = ''))
      print('This analysis had to be stopped because there is no more synapomorphies!')
      break
    }
    
    # roundNRes <- rownames(matTemp2[rowSums(matTemp2) < resulT, , drop = F])
    # roundNullRes <- rownames(matTemp2[rowSums(matTemp2) != 0, , drop = F])
    # roundInt <- intersect(roundNRes, roundNullRes)
    
    # ###############################################################
    # # após a primeira rodada e exclusão de caracteres com IC e IR iguais a 1...
    # lista_T <- list()
    # 
    # if(length(colnames(matTemp)) == 0){
    #   nomesCOL <- NULL
    #   stop('There is no more characters for the next iteration.')
    # }
    # 
    # nomesCOL <- colnames(matTemp) # species which are not synapomorphies
    # 
    # for(j in nomesCOL){
    #   discre <- matTemp[, as.character(j)]
    #   for(a in 1:length(discre)){
    #     ifelse(discre[a] == 1, discre[a] <- j, discre[a] <- 'absent')
    #   } # vector with the presence and absence of one species
    #   #This maps the character (i.e., taxon).
    #   # x11()
    #   # plotBranchbyTrait(chartree, tempMatrix[,which(ciVec == 1)])
    #   ## simulate single stochastic character map using empirical Bayes method
    #   discre <- as.factor(discre)  ### vetor de presença da espécie
    #   
    #   lista_T[[j]] <- names(discre)[which(discre == j)]
    # }
    # 
    # speciesNames <- names(lista_T)
    # 
    # qde <- 0
    # idx <- 1
    # for(i in speciesNames){
    #   idx <- idx + 1
    #   qde[idx] <- length(unlist(lista_T[speciesNames[idx - 1]]))
    # }
    # 
    # qdeCum <- cumsum(qde)
    # 
    # tabela <- matrix(NA, nrow = qdeCum[length(qde)], nc = 2)
    # 
    # indice <- 0
    # for(j in speciesNames){
    #   indice <- indice + 1
    #   tabela[(qdeCum[indice] + 1): qdeCum[indice + 1], 1] <- as.character(rep(j, qde[indice + 1]))
    #   tabela[(qdeCum[indice] + 1): qdeCum[indice + 1], 2] <- as.numeric(unlist(lista_T[[indice]]))
    # }
    # 
    # # tabela <- as.data.frame(tabela)
    # 
    # # data frame with the results...
    # frameTemp <- data.frame(spp = tabela[, 1], grid_n = tabela[, 2], row.names = 1:dim(tabela)[1])
    # # frameTemp
    # 
    # # matrix with species names and grid numbers:
    # resulPaeRaster <- matrix(as.matrix(frameTemp[,2]), nrow(frameTemp),
    #                          1, dimnames = list(frameTemp[,1], colnames(frameTemp)[2]))
    # 
    # 
    # ## the number of species supporting the grid number
    # speciesNumber <- data.frame(table(resulPaeRaster))
    # 
    # n_occur <- data.frame(table(unlist(lista_T))) #frequency of each grid
    # 
    # 
    # syn_grids <- unlist(lista_T)[unlist(lista_T) %in% n_occur$Var1[n_occur$Freq > 1]]
    # # syn_grids
    # if(length(syn_grids) < 2){
    #   print('This iteration has to be stopped because there is no more synapomorphies!')
    #   next
    # }
    ################################################################
    ###########################
    ## para ver o mapeamento por localidade das espécies que formam
    # o(s) traço(s) genberalizado(s):
    # após a primeira rodada e exclusão de caracteres com IC e IR iguais a 1...
    lista_TT <- list()
    
    # nomesCOLsim # species which are synapomorphies (first stage)
    # roundRes # número de grid que são os traços generalizados
    
    # common species which are synapomorphies
    nLi <- dim(matTemp2)[1]
    roundResCol1 <- colnames(matTemp2[, colSums(matTemp2) <= nLi , drop = F])
    roundResCol2 <- colnames(matTemp2[, colSums(matTemp2) > 1 , drop = F])
    roundResCol <- intersect(roundResCol1, roundResCol2)
    
    if(is.null(roundResCol) && contagem == 1){
      plot(cropped_map, add = TRUE)
      stop('This analysis has no indicative of generalized tracks!')
    } else if(is.null(roundResCol) && contagem > 1){
      print(paste0(c('The iteration number ', contagem, ' has finished without results.'), collapse = ''))
      print('This analysis had to be stopped because there is no more synapomorphies!')
      break
    }
    
    if(length(roundResCol) < 2 && contagem == 1){
      plot(cropped_map, add = TRUE)
      stop('This analysis has no indicative of generalized tracks!')
    } else if(length(roundResCol) < 2 && contagem > 1){
      print(paste0(c('The iteration number ', contagem, ' has finished without results.'), collapse = ''))
      print('This analysis had to be stopped because there is no more synapomorphies!')
      break
    }
    
    quantum <- 1
    for(j in 1:length(roundResCol)){
      k <- roundResCol[quantum]
      discre <- matTemp2[, as.character(k)]
      for(a in 1:length(discre)){
        ifelse(discre[a] == 1, discre[a] <- k, discre[a] <- 'absent')
      } # vector with the presence and absence of one species
      #This maps the character (i.e., taxon).
      # x11()
      # plotBranchbyTrait(chartree, tempMatrix[,which(ciVec == 1)])
      ## simulate single stochastic character map using empirical Bayes method
      
      if(all(discre == 'absent')){
        quantum <- quantum + 1
        print(paste0(k, 'is outside of map! Please check it!'))
        next
      }
      
      discre <- as.factor(discre)  ### vetor de presença da espécie
      
      lista_TT[[k]] <- names(discre)[which(discre == k)]
      
      quantum <- quantum + 1
    }
    
    # adic <- dim(matTemp2)[2] - length(which(colnames(matTemp2) == roundResCol))
    subtrair <- NULL
    for(i in colnames(matTemp2)){
      subtrair[i] <- which(colnames(preabsMat) == i)
    }
    
    # posic <- 0
    # for(pos in 1:length(unique(rownames(resulPaeRaster)))){
    #   posic[pos] <- which(colnames(preabsMat) == unique(rownames(resulPaeRaster))[pos])
    # }
    # for(j in roundResCol){
    #   subtrair <- subtrair[-which(names(subtrair) == j)]
    # }
    # 
    # for(i in 1:length(subtrair)){
    #   posic[pos + i] <- which(colnames(preabsMat) == names(subtrair[i]))
    # }
    
    # changing the matrix for the next iteration:	
    tempMatrix <- subset(tempMatrix, select=-subtrair, site.pattern = FALSE)
    
    # matTemp <- as.data.frame(matTemp)
    # tamanhoMat <- dim(matTemp)[2]
    # matTemp <- cbind(matTemp, preabsMat[, subtrair])
    
    # colnames(matTemp) <- c(colnames(matTemp[, 1:tamanhoMat]), names(subtrair)) 
    
    speciesNamesT <- names(lista_TT)
    
    qde <- 0
    idx <- 1
    for(i in speciesNamesT){
      idx <- idx + 1
      qde[idx] <- length(unlist(lista_TT[speciesNamesT[idx - 1]]))
    }
    
    qdeCum <- cumsum(qde)
    
    tabelaT <- matrix(NA, nrow = qdeCum[length(qde)], nc = 2)
    
    indice <- 0
    for(j in speciesNamesT){
      indice <- indice + 1
      tabelaT[(qdeCum[indice] + 1): qdeCum[indice + 1], 1] <- as.character(rep(j, qde[indice + 1]))
      tabelaT[(qdeCum[indice] + 1): qdeCum[indice + 1], 2] <- as.numeric(unlist(lista_TT[[indice]]))
    }
    
    # tabela <- as.data.frame(tabela)
    
    if(dim(tabelaT)[1] == 0){
      stop('No more iterations are needed.')
    }
    # data frame with the results...
    frameTempT <- data.frame(spp = tabelaT[, 1], grid_n = tabelaT[, 2], row.names = 1:dim(tabelaT)[1])
    # frameTemp
    
    # matrix with species names and grid numbers:
    resulPaeRasterT <- matrix(as.matrix(frameTempT[,2]), nrow(frameTempT),
                             1, dimnames = list(frameTempT[,1], colnames(frameTempT)[2]))
    
    
    ## the number of species supporting the grid number
    speciesNumberT <- data.frame(table(resulPaeRasterT))
    
    n_occurT <- data.frame(table(unlist(lista_TT))) #frequency of each grid
    
    similarNames <- list()
    for(volta in 1:dim(speciesNumberT)[1]){
      similarNames[[volta]] <- resulPaeRasterT[which(resulPaeRasterT[,1] == frameTempT$grid_n[volta]), ]
      names(similarNames[[volta]]) <- subset(frameTempT, grid_n == similarNames[[volta]][1])[,1]
    }
    
    ## grid numbers that are generalized tracks in this iteration
    gridIt <- unique(as.numeric(unlist(similarNames)))
    # gridIt
    
    
    if(length(unique(rownames(resulPaeRasterT))) == 0){
      stop('This iteration had to be stopped because there is no more synapomorphies!')
    }
    
    if(length(unique(frameTempT$spp)) > 1){
      # set the cells associated with the shapfile to the specified value
      r[r > 0] <- NA
      
      for(j in 1:dim(speciesNumberT)[1]){
        if(speciesNumberT$Freq[j] > 1){
          gTrack <- as.numeric(as.character(speciesNumberT[j,1]))
          values(r)[gTrack] <- 1
        }
      }
      
      if(sobrepo == TRUE){
        plot(r, axes = FALSE, legend = FALSE, add = TRUE, col = colores[contagem],
             alpha = 0.60)
      } else if(sobrepo == FALSE){
        plot(r, axes = FALSE, legend = FALSE, add = TRUE, col = colores[1],
             alpha = 0.50)
      }
    
      
      # producing a raster:
      #convert the raster to points for plotting the number of a grid
      map.r <- as.data.frame(rasterToPoints(r))
      pontosRaster <- rasterize(cbind(map.r$x, map.r$y), r, field = 1) # raster com as presenças
      proj4string(pontosRaster) <- CRS("+proj=longlat +datum=WGS84")
      writeRaster(pontosRaster, paste0(c('out/generalizedTrack_', contagem, '.tif'), collapse = ''),
                      overwrite = TRUE)
      
      
      ## nota: gridPolygon@data has all polygons established with a specific resolution
      
      if(gridView == TRUE){
        plot(cropped_map, add = TRUE)
        
        map.r$gridNumber <- which(pontosRaster@data@values == 1)
        
        if(labelGrid == TRUE){
          text(map.r[,c(1, 2)], labels = map.r$gridNumber, cex = 2.0, col = 'black', font = 2)
        }
      }
      
      if(legendSpecies == TRUE){
        print('Please, check the legendSpeciesA.txt file in the /out directory...')
        #local <- locator()
        logfile <- "out/legendSpeciesA.txt"
        cat(c("grid_number", "species", "\n"), file = logfile, sep="\t")
        for(i in as.numeric(as.character(speciesNumberT$resulPaeRasterT))){
          cat(c(paste0(c('grid_', i), collapse = ''), rownames(subset(resulPaeRasterT,
                               resulPaeRasterT[,'grid_n'] == i)), '\n'), sep = '\t', file = logfile, append = TRUE)
        }
      }
    }
    # 
    # ##################################
    # ## Mudando matTemp:
    ######################################################
    # interT2 <- intersect(interT1, inter3)
    matTemp <- preabsMat[, -subtrair] # temporary matrix with homoplastic synapomorphies
    matTemp <- as.matrix(matTemp)
    # colnames(matTemp) <- colnames(preabsMat)
    #######################################################
    
    
    nomesCOL <- colnames(matTemp) # species which are not synapomorphies

    quantum <- 1
    for(j in 1:length(nomesCOL)){
      k <- nomesCOL[quantum]
      discre <- matTemp[, as.character(k)]
      for(a in 1:length(discre)){
        ifelse(discre[a] == 1, discre[a] <- k, discre[a] <- 'absent')
      } # vector with the presence and absence of one species
      #This maps the character (i.e., taxon).
      # x11()
      # plotBranchbyTrait(chartree, tempMatrix[,which(ciVec == 1)])
      ## simulate single stochastic character map using empirical Bayes method
      
      if(all(discre == 'absent')){
        quantum <- quantum + 1
        print(paste0(k, 'is outside of map! Please check it!'))
        break
      }
      
      discre <- as.factor(discre)  ### vetor de presença da espécie
      
      lista_T[[k]] <- names(discre)[which(discre == k)]
      
      quantum <- quantum + 1
      # print(length(unlist(lista_T[[k]])))
    }
    
    # unlist(lista_T)[unlist(lista_T) %in% n_occur$Var1[n_occur$Freq > 1]]

    speciesNames <- names(lista_T)

    qde <- 0
    idx <- 1
    for(i in speciesNames){
      idx <- idx + 1
      qde[idx] <- length(unlist(lista_T[speciesNames[idx - 1]]))
    }

    qdeCum <- cumsum(qde)

    tabela <- matrix(NA, nrow = qdeCum[length(qde)], nc = 2)

    indice <- 0
    for(j in speciesNames){
      indice <- indice + 1
      tabela[(qdeCum[indice] + 1): qdeCum[indice + 1], 1] <- as.character(rep(j, qde[indice + 1]))
      tabela[(qdeCum[indice] + 1): qdeCum[indice + 1], 2] <- as.numeric(unlist(lista_T[[indice]]))
    }

    # tabela <- as.data.frame(tabela)

    if(dim(tabela)[1] == 0){
      print('No more iterations are needed.')
      break
    }

    # data frame with the results...
    frameTemp <- data.frame(spp = tabela[, 1], grid_n = tabela[, 2], row.names = 1:dim(tabela)[1])
    # frameTemp

    # matrix with species names and grid numbers:
    resulPaeRaster <- matrix(as.matrix(frameTemp[,2]), nrow(frameTemp),
                             1, dimnames = list(frameTemp[,1], colnames(frameTemp)[2]))


    ## the number of species supporting the grid number
    speciesNumber <- data.frame(table(resulPaeRaster))

    n_occur <- data.frame(table(unlist(lista_T))) #frequency of each grid


    syn_grids <- unlist(lista_T)[unlist(lista_T) %in% n_occur$Var1[n_occur$Freq > 1]]
    # syn_grids
    if(length(syn_grids) < 2){
      print('This iteration has to be stopped because there is no more synapomorphies!')
      next
    }
    
    lista[[contagem]] <- frameTemp  # species which are not synapomorphies
    listaR[[contagem]] <- frameTempT # species which are synapomorphies
    
    
    # if(length(roundResCol) < 2){
    #   
    #   next('This iteration has only one synapomorphy...')
    # }
    
    preabsMat <- matTemp
    
    print(paste0(c('The iteration number ', contagem, ' has finished.'), collapse = ''))

    contagem <- contagem + 1
    conta[contagem] <- contagem
  } # close while looping
  
  
  if(length(colnames(matTemp)) == 0 && length(conta) == 1){
    xis <- 1
    plot(cropped_map, add = TRUE)
    stop('This analysis has no indicative of generalized tracks!')
  } else if(length(colnames(matTemp)) == 0 && length(conta) > 1){
    xis <- seq(1:length(listaR))
  }
  
  if(is.null(nomesCOL) || dim(as.data.frame(tempMatrix))[1] == 0){
    xis <- seq(1:length(listaR))
  }
  
  
  if(labelGrid == TRUE){
    if(sobrepo == TRUE){
      print('Bold numbers refer to PAE analysis using generalized tracks!')
    }
    text(map.r[,c(1, 2)], labels = map.r$gridNumber, cex = 2.0, col = 'black', font = 2)
  }
  
  conta <- NULL
  if(sobrepo == FALSE){
    # for(i in 1:contagem){
    #   if(length(unique(listaR[[i]]$spp)) == 1){
    #     conta[i] <- 'singleton'
    #   } else if(length(unique(listaR[[i]]$spp)) > 1){
    #     conta[i] <- 'non_singleton'
    #   }
    # }
    # xis <- as.numeric(which(conta == 'non_singleton'))
    xis <- seq(1:length(listaR))
    legend(x = 'topright', legend = xis, pch = 15, col = colores[xis],
           title = 'Adding generalized track of the iterations...', title.col = 'red', pt.cex = 1.5, cex = 0.8)
  } else if(sobrepo == TRUE){
    # for(i in 1:contagem){
    #   if(length(unique(listaR[[i]]$spp)) == 1){
    #     conta[i] <- 'singleton'
    #   } else if(length(unique(listaR[[i]]$spp)) > 1){
    #     conta[i] <- 'non_singleton'
    #   }
    # }
    # xis <- as.numeric(which(conta == 'non_singleton'))
    xis <- seq(1:length(listaR))
    legend(x = 'topright', legend = xis, pch = 15, col = colores[xis],
           title = 'Adding generalized track of the iterations...', title.col = 'red', pt.cex = 1.5, cex = 0.8)
  }
  
  if(contagem == 1){
    print('There is no result!')
  } else if(contagem != 1){
    return(list(homoplastic_species = lista,
                nonHomoplastic_species = listaR))
  }
} # close if... else... condition