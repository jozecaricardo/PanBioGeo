#### PAE-PCE analysis - VERSION WITH CORRECT STOPPING CRITERIA ####
# preabsMat = the presence-absence matrix from terminal_node function
# rou = number of iterations in the parsimony with N runs
pae_pce <- function(preabsMat, shapeFile, resolut, N = NULL,
                    gridView = FALSE, labelGrid = FALSE, legendSpecies = TRUE, sobrepo = FALSE,
                    xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
                    min_sobreposicao = 1.0){
  
  library(phangorn)
  library(raster)
  library(viridis)
  library(TreeSearch) # TBR search
  
  ## HELPER FUNCTION FOR ELEGANT TERMINATIONS ##
  finalizar_analise_elegante <- function(reason, details = NULL, iteration = 1, plot_grid = TRUE) {
    # Function to gracefully finalize the analysis with informative reports
    
    print("")
    print("")
    print(paste(rep("=", 60), collapse=""))
    print("FINALIZED PAE-PCE ANALYSIS")
    print(paste(rep("=", 60), collapse=""))
    
    print(paste0("Stopping iteration: ", iteration))
    print(paste0("Reason for termination: ", reason))
    
    if(!is.null(details)) {
      print("Details:")
      for(detail in details) {
        print(paste0("  • ", detail))
      }
    }
    
    print("")
    print("Scientific interpretation:")
    if(grepl("synapomorphy", reason, ignore.case = TRUE)) {
      print("• This is a normal completion of the PAE-PCE analysis")
      print("• Indicates that there are no more informative characters for analysis")
      print("• All non-homoplastic synapomorphies were identified")
    } else if(grepl("overlap", reason, ignore.case = TRUE)) {
      print("• The analyzed species do not have sufficiently overlapping distributions")
      print("• Consider adjusting the resolution of the grids or the overlap criteria")
      print("• It may indicate that the species have naturally distinct distributions")
    } else {
      print("• This is a normal termination based on the analysis criteria")
      print("• The data do not meet the minimum requirements for generalized traits")
    }
    
    print("")
    print("Recommendations:")
    print("• Check the quality of the input data")
    print("• Consider adjusting resolution parameters")
    print("• Assess whether the data are suitable for PAE-PCE analysis")
    
    print(paste(rep("=", 60), collapse=""))
    
    # Plot grid if requested
    if(plot_grid && exists("cropped_map")) {
      plot(cropped_map, add = TRUE, border = "gray", lwd = 0.5)
    }
    
    # Return structured information instead of throwing error
    return(list(
      status = "finished_without_results",
      reason = reason,
      details = details,
      iteration_parada = iteration,
      timestamp = Sys.time(),
      message = "PAE-PCE analysis completed - please refer to the report for details"
    ))
  }
  
  ## HELPER FUNCTION TO CHECK GRID OVERLAP - CRITERION 1.0 ##
  verificar_sobreposicao_grids <- function(matTemp2, min_sobreposicao = 1.0) {
    # Identifies groups of species with PERFECT overlap (1.0 = identical distributions)
    
    print("=== CHECKING GROUPS WITH PERFECT OVERLAP ===")
    
    # Get list of grids for each species
    especies_grids <- list()
    for(especie in colnames(matTemp2)) {
      grids_especie <- rownames(matTemp2)[matTemp2[, especie] == 1]
      especies_grids[[especie]] <- grids_especie
      print(paste0(especie, ": grids ", paste(grids_especie, collapse = ", ")))
    }
    
    especies <- names(especies_grids)
    n_especies <- length(especies)
    
    if(n_especies < 2) {
      print("Only 1 species - automatically accepted")
      return(list(especies_validas = especies, sobreposicao_ok = TRUE))
    }
    
    # Calculate overlap matrix among all species
    sobreposicoes <- matrix(0, nrow = n_especies, ncol = n_especies, 
                            dimnames = list(especies, especies))
    
    for(i in 1:n_especies) {
      sobreposicoes[i, i] <- 1  # Self-overlap = 1
    }
    
    for(i in 1:(n_especies-1)) {
      for(j in (i+1):n_especies) {
        esp1 <- especies[i]
        esp2 <- especies[j]
        grids1 <- especies_grids[[esp1]]
        grids2 <- especies_grids[[esp2]]
        
        # Calculate overlap (Jaccard index)
        intersecao <- length(intersect(grids1, grids2))
        uniao <- length(union(grids1, grids2))
        sobreposicao <- if(uniao > 0) intersecao / uniao else 0
        
        sobreposicoes[i, j] <- sobreposicao
        sobreposicoes[j, i] <- sobreposicao
        
        print(paste0("Overlap ", esp1, " vs ", esp2, ": ", round(sobreposicao, 3)))
        
        # STRICT CRITERION: Only overlap = 1.0 is accepted
        if(sobreposicao == 1.0) {
          print(paste0("  ✓ PERFECT OVERLAP (= 1.0) - Identical distributions"))
        } else {
          print(paste0("  ✗ Imperfect overlap (", round(sobreposicao, 3), " ≠ 1.0) - Different distributions"))
        }
        print(paste0("  - Common grids: ", paste(intersect(grids1, grids2), collapse = ", ")))
        if(sobreposicao < 1.0) {
          print(paste0("  - Unique grids of ", esp1, ": ", paste(setdiff(grids1, grids2), collapse = ", ")))
          print(paste0("  - Unique grids of ", esp2, ": ", paste(setdiff(grids2, grids1), collapse = ", ")))
        }
      }
    }
    
    # LOGIC: Find groups of species with PERFECT overlap (1.0)
    print("")
    print("=== IDENTIFYING GROUPS WITH IDENTICAL DISTRIBUTIONS ===")
    
    # Create connection graph (species connected ONLY if overlap = 1.0)
    conexoes <- sobreposicoes == 1.0
    
    # Find connected components (species groups)
    visitado <- rep(FALSE, n_especies)
    grupos <- list()
    grupo_id <- 1
    
    for(i in 1:n_especies) {
      if(!visitado[i]) {
        # Depth-first search to find all connected species
        grupo_atual <- c()
        pilha <- c(i)
        
        while(length(pilha) > 0) {
          atual <- pilha[length(pilha)]
          pilha <- pilha[-length(pilha)]
          
          if(!visitado[atual]) {
            visitado[atual] <- TRUE
            grupo_atual <- c(grupo_atual, atual)
            
            # Add unvisited neighbors to stack
            vizinhos <- which(conexoes[atual, ] & !visitado)
            pilha <- c(pilha, vizinhos)
          }
        }
        
        if(length(grupo_atual) >= 2) {
          grupos[[grupo_id]] <- especies[grupo_atual]
          print(paste0("Group ", grupo_id, " (", length(grupo_atual), " species with identical distributions): ", 
                       paste(especies[grupo_atual], collapse = ", ")))
          grupo_id <- grupo_id + 1
        } else {
          print(paste0("Isolated species (unique distribution): ", especies[grupo_atual]))
        }
      }
    }
    
    # Determine valid species (those belonging to groups ≥2 species with overlap = 1.0)
    especies_validas <- c()
    for(grupo in grupos) {
      especies_validas <- c(especies_validas, grupo)
    }
    
    if(length(especies_validas) >= 2) {
      print(paste0("", "✓ FOUND ", length(grupos), " GROUP(S) WITH IDENTICAL DISTRIBUTIONS"))
      print(paste0("Accepted species (identical distributions): ", paste(especies_validas, collapse = ", ")))
      return(list(especies_validas = especies_validas, sobreposicao_ok = TRUE,
                  sobreposicoes = sobreposicoes, grupos = grupos))
    } else {
      print("")
      print("✗ NO VALID GROUP FOUND (no group with ≥2 species of identical distributions)")
      return(list(especies_validas = c(), sobreposicao_ok = FALSE,
                  sobreposicoes = sobreposicoes, grupos = list()))
    }
  }
  
  # Convert the character matrix to a phyDat object 
  # in order to infer a tree. By assigning the matrix 
  # type as 'user' we can specify the components of the 
  # matrix (in this case, binary 1s and 0s).
  
  tempMatrix <- as.phyDat(preabsMat, type = 'USER', levels = c(0, 1))
  
  colu <- ncol(preabsMat)
  
  # Creating a restriction:
  ciVec <- rep(1.0, colu)
  riVec <- rep(1.0, colu)
  
  
  lista <- list()
  listaR <- list()
  
  print('Please, remember that the resulting rasters will be kept in your setted directory!')
  
  if(sobrepo == TRUE){
    print('Please do not forgive that grids are in different colors!')
    if(is.null(N)){
      # ELEGANT STOP: Mandatory parameter not provided
      return(finalizar_analise_elegante(
        reason = "Parameter N not specified",
        details = c("Parameter N (number of iterations) is mandatory when sobrepo=TRUE",
                    "Specify N=10 or another suitable value"),
        iteration = 0,
        plot_grid = FALSE
      ))
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
  
  
  # CORRECTION: Creating cropped_map with perfect cut at shapefile boundaries
  # First create the base grid
  grid_base <- raster(extent(shapeFile), resolution = resolut, crs = CRS("+proj=longlat +datum=WGS84"))
  grid_base <- rasterize(shapeFile, grid_base, field = 1)
  
  # Convert to polygons only cells that are within the shapefile
  grid_base[is.na(grid_base)] <- 0
  grid_base[grid_base == 0] <- NA
  cropped_map <- rasterToPolygons(grid_base, dissolve = FALSE)
  
  # Ensure it has the same projection
  crs(cropped_map) <- "+proj=longlat +datum=WGS84"
  
  # Producing a raster of the shapefile
  mask.raster <- raster(extent(shapeFile), resolution = resolut,
                        crs = CRS("+proj=longlat +datum=WGS84"))
  r <- rasterize(shapeFile, mask.raster)
  proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  r <- merge(r, mask.raster)
  
  
  if(sobrepo == FALSE){
    colores <- viridis(n = length(r[!is.na(r)]), option = 'A')
  } else if(sobrepo == TRUE){
    colores <- viridis(n = N, option = 'A')
  }
  
  # Plotting the shapefile
  if(is.null(xmin)){
    plot(shapeFile, axes = TRUE, las = 1)
  } else if(!is.null(xmin)){
    plot(shapeFile, axes = TRUE, las = 1, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  }
  
  
  #---------------------------------- #### PAE-PCE ### -------------------------------------#
  contagem <- 1
  conta <- 0
  
  # MAIN CORRECTION: Variable to control when to stop the loop
  deve_continuar <- TRUE
  
  # CORRECTION: More robust while condition
  while(deve_continuar && length(which((ciVec == 1) == TRUE)) > 1 && dim(preabsMat)[2] > 1){
    
    print(paste0("=== STARTING ITERATION ", contagem, " ==="))
    print(paste0("Remaining species in matrix: ", dim(preabsMat)[2]))
    print(paste0("Grids in matrix: ", dim(preabsMat)[1]))
    
    # We can get a random starting tree for a parsimony search using the
    # random.addition function:
    
    if(dim(as.data.frame(tempMatrix))[1] < 2){
      # ELEGANT STOP: Insufficient data
      return(finalizar_analise_elegante(
        reason = "Insufficient data for phylogenetic analysis",
        details = c("Less than 2 terminals in the data matrix",
                    "Parsimony analysis requires at least 2 terminals"),
        iteration = contagem
      ))
    }
    
    ### TBR + RATCHET searches
    dm <- dist.hamming(tempMatrix)
    ra.tre <- NJ(dm)
    
    # Use the pratchet function (parsimony ratchet) to find the most
    # parsimonious tree. Specifying k means the algorithm will search   
    # through k possible trees to find the most parsimonious solution.
    treeIt1 <- list()
    
    for(rou in 1:N){
      treeIt1[[rou]] <- optim.parsimony(tree = ra.tre, data = tempMatrix,
                                        method = 'fitch', rearrangements = 'SPR')
      treeIt1[[rou]] <- pratchet(tempMatrix, start = treeIt1[[rou]], maxit=1000,
                                 minit=100, k=10, trace=1)
      
      treeIt1[[rou]] <- consensus(treeIt1[[rou]])
      
      ra.tre <- treeIt1[[rou]]
    }
    
    if(is.null(treeIt1)){
      # ELEGANT STOP: Tree construction failure
      return(finalizar_analise_elegante(
        reason = "Failure in phylogenetic topology construction",
        details = c("Search algorithm could not generate valid trees",
                    "May indicate problems in the input data"),
        iteration = contagem
      ))
    }
    
    
    
    treeT <- consensus(treeIt1) # 'majority rule': change p to 0.5
    
    # Root the tree by the designated outgroup 
    # (write the species name as it appears in the tree, 
    # but add an underscore to fill the space).
    strictTree <- root(treeT, outgroup = "ROOT", resolve.root = TRUE)
    
    # CI of each character (i.e., taxon):
    ciVec <- CI(strictTree, data = tempMatrix, sitewise = TRUE)
    
    # RI of each character (i.e., taxon):
    riVec <- RI(strictTree, data = tempMatrix, sitewise = TRUE)
    riVec[which(is.na(riVec))] <- 0
    
    print(paste0("CI calculated for ", length(ciVec), " characters"))
    print(paste0("RI calculated for ", length(riVec), " characters"))
    print(paste0("Characters with CI=1: ", length(which(ciVec == 1))))
    print(paste0("Characters with RI=1: ", length(which(riVec == 1))))
    
    if(all(is.na(riVec) == TRUE)){
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break # remaining body of loop will not be executed!
    }
    
    if(any(riVec == 1.0) == FALSE){
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break # remaining body of loop will not be executed!
    }
    
    
    ## First load packages
    require(phytools)
    tree_temp <- strictTree
    if(is.null(tree_temp$edge.length) == TRUE){
      print('...creating the branch lengths of a tree equal to one...')
      tree_temp <- compute.brlen(tree_temp, 1)
    }
    
    # Discrete mapping
    inter1 <- which(ciVec == 1)
    inter2 <- which(riVec == 1)
    interT1 <- intersect(inter1, inter2)
    
    print(paste0("Characters with CI=1 and RI=1 (non-homoplastic synapomorphies): ", length(interT1)))
    
    if(length(interT1) == 0 && contagem == 1){
      # ELEGANT STOP: No non-homoplastic synapomorphy in first iteration
      return(finalizar_analise_elegante(
        reason = "No non-homoplastic synapomorphy found",
        details = c("First iteration did not identify characters with CI=1 and RI=1",
                    "Data may not be suitable for PAE-PCE analysis",
                    "Consider checking the quality of the input data"),
        iteration = contagem
      ))
    }
    
    if(length(interT1) < 2){
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break # remaining body of loop will not be executed!
    }
    
    # Things that had IC and IR equal to 1...
    lista_T <- list()
    matTemp2 <- preabsMat[, interT1] # temporary matrix with non-homoplastic synapomorphies
    matTemp2 <- as.matrix(matTemp2)
    
    nomesCOLsim <- colnames(matTemp2) # species which are synapomorphies
    print(paste0("Synapomorphic species identified: ", paste(nomesCOLsim, collapse = ", ")))
    
    if(length(nomesCOLsim) < 2){
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break # remaining body of loop will not be executed!
    }
    
    resulT <- dim(matTemp2)[2]
    
    # NEW VERIFICATION: Strict overlap of distributions
    print(paste0('Iteration ', contagem, ': Identifying groups of species with identical distributions...'))
    
    verificacao <- verificar_sobreposicao_grids(matTemp2, min_sobreposicao)
    
    # CRITICAL CORRECTION: Stop iteration if no perfect overlap found
    if(!verificacao$sobreposicao_ok || length(verificacao$especies_validas) < 2) {
      # No valid overlap group found
      if(contagem == 1) {
        # First iteration - elegant stop
        return(finalizar_analise_elegante(
          reason = "No valid group of species with overlapping distributions",
          details = c("First iteration found no species with identical distributions",
                      paste0("Species analyzed: ", paste(colnames(matTemp2), collapse = ", ")),
                      "Consider adjusting grid resolution or overlap criteria"),
          iteration = contagem
        ))
      } else {
        # Subsequent iterations - STOP HERE, DO NOT CONTINUE
        print("")
        print(paste(rep("=", 60), collapse=""))
        print(paste0("ITERATION ", contagem, ": NO PERFECT OVERLAP FOUND"))
        print(paste(rep("=", 60), collapse=""))
        print("No group with identical distributions (overlap = 1.0) was found.")
        print("Analysis MUST stop here - no further iterations allowed.")
        print(paste0("Species analyzed: ", paste(colnames(matTemp2), collapse = ", ")))
        print(paste(rep("=", 60), collapse=""))
        deve_continuar <- FALSE
        break  # STOP THE LOOP IMMEDIATELY
      }
    }
    
    # Filter only species with identical distributions
    especies_validas <- verificacao$especies_validas
    matTemp2 <- matTemp2[, especies_validas, drop = FALSE]
    
    print(paste0('Accepted species (identical distributions): ', paste(especies_validas, collapse = ", ")))
    
    # Continue with normal processing
    
    # Common grids (synapomorphies)
    roundRes1 <- rownames(matTemp2[rowSums(matTemp2) <= resulT, , drop = F])
    roundRes2 <- rownames(matTemp2[rowSums(matTemp2) > 1, , drop = F])
    roundRes <- intersect(roundRes1, roundRes2)
    
    print(paste0("Common grids identified: ", length(roundRes)))
    if(length(roundRes) > 0) {
      print(paste0("Grids: ", paste(roundRes, collapse = ", ")))
    }
    
    if(is.null(roundRes) && contagem == 1){
      # ELEGANT STOP: No common grid in first iteration
      return(finalizar_analise_elegante(
        reason = "No common grid identified as synapomorphy",
        details = c("First iteration found no shared grids among species",
                    "Species may have very dispersed distributions",
                    "Consider adjusting grid resolution"),
        iteration = contagem
      ))
    } else if(is.null(roundRes) && contagem > 1){
      print(paste0(c('The iteration number ', contagem, ' has finished without results.'), collapse = ''))
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break
    }
    
    ################################################################
    ###########################
    ## To see the mapping by locality of the species that form
    # the generalized track(s):
    # After the first round and exclusion of characters with IC and IR equal to 1...
    lista_TT <- list()
    
    # Common species which are synapomorphies
    nLi <- dim(matTemp2)[1]
    roundResCol1 <- colnames(matTemp2[, colSums(matTemp2) <= nLi , drop = F])
    roundResCol2 <- colnames(matTemp2[, colSums(matTemp2) > 1 , drop = F])
    roundResCol <- intersect(roundResCol1, roundResCol2)
    
    print(paste0("Common synapomorphic species identified: ", length(roundResCol)))
    if(length(roundResCol) > 0) {
      print(paste0("Species: ", paste(roundResCol, collapse = ", ")))
    }
    
    if(is.null(roundResCol) && contagem == 1){
      # ELEGANT STOP: No common synapomorphic species in first iteration
      return(finalizar_analise_elegante(
        reason = "No common synapomorphic species identified",
        details = c("First iteration found no shared species among grids",
                    "May indicate that species have very specific distributions",
                    "Consider checking species distribution data"),
        iteration = contagem
      ))
    } else if(is.null(roundResCol) && contagem > 1){
      print(paste0(c('The iteration number ', contagem, ' has finished without results.'), collapse = ''))
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break
    }
    
    if(length(roundResCol) < 2 && contagem == 1){
      # ELEGANT STOP: Less than 2 synapomorphic species in first iteration
      return(finalizar_analise_elegante(
        reason = "Insufficient synapomorphic species for analysis",
        details = c("First iteration found less than 2 synapomorphic species",
                    paste0("Species found: ", length(roundResCol)),
                    "PAE-PCE analysis requires at least 2 synapomorphic species"),
        iteration = contagem
      ))
    } else if(length(roundResCol) < 2 && contagem > 1){
      print(paste0(c('The iteration number ', contagem, ' has finished without results.'), collapse = ''))
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break
    }
    
    quantum <- 1
    for(j in 1:length(roundResCol)){
      k <- roundResCol[quantum]
      discre <- matTemp2[, as.character(k)]
      for(a in 1:length(discre)){
        ifelse(discre[a] == 1, discre[a] <- k, discre[a] <- 'absent')
      } # vector with the presence and absence of one species
      
      if(all(discre == 'absent')){
        quantum <- quantum + 1
        print(paste0(k, 'is outside of map! Please check it!'))
        next
      }
      
      discre <- as.factor(discre)  ### species presence vector
      
      lista_TT[[k]] <- names(discre)[which(discre == k)]
      
      quantum <- quantum + 1
    }
    
    subtrair <- NULL
    for(i in colnames(matTemp2)){
      subtrair[i] <- which(colnames(preabsMat) == i)
    }
    
    print(paste0("Removing ", length(subtrair), " species from matrix for next iteration"))
    print(paste0("Species to remove: ", paste(names(subtrair), collapse = ", ")))
    
    # Changing the matrix for the next iteration:	
    tempMatrix <- subset(tempMatrix, select=-subtrair, site.pattern = FALSE)
    
    print(paste0("New tempMatrix: ", dim(as.data.frame(tempMatrix))[2], " species"))
    
    # CRUCIAL VERIFICATION: If there are not enough species, stop
    if(dim(as.data.frame(tempMatrix))[2] < 2) {
      print("Less than 2 species remaining in matrix - stopping analysis")
      deve_continuar <- FALSE
      break
    }
    
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
    
    if(dim(tabelaT)[1] == 0){
      # ELEGANT STOP: Empty results table
      return(finalizar_analise_elegante(
        reason = "No results generated in the iteration",
        details = c("Species and grids table resulted empty",
                    "May indicate natural endpoint of PAE-PCE analysis"),
        iteration = contagem
      ))
    }
    # Data frame with the results...
    frameTempT <- data.frame(spp = tabelaT[, 1], grid_n = tabelaT[, 2], row.names = 1:dim(tabelaT)[1])
    
    # Matrix with species names and grid numbers:
    resulPaeRasterT <- matrix(as.matrix(frameTempT[,2]), nrow(frameTempT),
                              1, dimnames = list(frameTempT[,1], colnames(frameTempT)[2]))
    
    
    ## The number of species supporting the grid number
    speciesNumberT <- data.frame(table(resulPaeRasterT))
    
    n_occurT <- data.frame(table(unlist(lista_TT))) # frequency of each grid
    
    similarNames <- list()
    for(volta in 1:dim(speciesNumberT)[1]){
      similarNames[[volta]] <- resulPaeRasterT[which(resulPaeRasterT[,1] == frameTempT$grid_n[volta]), ]
      names(similarNames[[volta]]) <- subset(frameTempT, grid_n == similarNames[[volta]][1])[,1]
    }
    
    ## Grid numbers that are generalized tracks in this iteration
    gridIt <- unique(as.numeric(unlist(similarNames)))
    
    
    if(length(unique(rownames(resulPaeRasterT))) == 0){
      # ELEGANT STOP: No unique species in result
      return(finalizar_analise_elegante(
        reason = "No unique species identified in result",
        details = c("Iteration result contains no distinct species",
                    "May indicate end of analysis or data problems"),
        iteration = contagem
      ))
    }
    
    if(length(unique(frameTempT$spp)) > 1){
      # Set the cells associated with the shapfile to the specified value
      r[r > 0] <- NA
      
      for(j in 1:dim(speciesNumberT)[1]){
        if(speciesNumberT$Freq[j] > 1){
          gTrack <- as.numeric(as.character(speciesNumberT[j,1]))
          values(r)[gTrack] <- 1
        }
      }
      
      # CORRECTION: Cut plotted cells at shapefile boundaries
      r_masked <- mask(r, shapeFile)
      
      if(sobrepo == TRUE){
        plot(r_masked, axes = FALSE, legend = FALSE, add = TRUE, col = colores[contagem],
             alpha = 0.60)
      } else if(sobrepo == FALSE){
        plot(r_masked, axes = FALSE, legend = FALSE, add = TRUE, col = colores[1],
             alpha = 0.50)
      }
      
      
      # Producing a raster:
      # Convert the raster to points for plotting the number of a grid
      map.r <- as.data.frame(rasterToPoints(r_masked))
      pontosRaster <- rasterize(cbind(map.r$x, map.r$y), r_masked, field = 1) # raster with presences
      proj4string(pontosRaster) <- CRS("+proj=longlat +datum=WGS84")
      writeRaster(pontosRaster, paste0(c('out/generalizedTrack_', contagem, '.tif'), collapse = ''),
                  overwrite = TRUE)
      
      
      ## Note: gridPolygon@data has all polygons established with a specific resolution
      
      if(gridView == TRUE){
        plot(cropped_map, add = TRUE, border = "gray", lwd = 0.5)
        
        map.r$gridNumber <- which(pontosRaster@data@values == 1)
        
        if(labelGrid == TRUE){
          text(map.r[,c(1, 2)], labels = map.r$gridNumber, cex = 0.8, col = 'black', font = 2)
        }
      }
      
      if(legendSpecies == TRUE){
        print('Please, check the legendSpeciesA.txt file in the out/ directory...')
        if (!file.exists('out/'))
          dir.create('out/')
        logfile <- "out/legendSpeciesA.txt"
        cat(c("grid_number", "species", "\n"), file = logfile, sep="\t")
        for(i in as.numeric(as.character(speciesNumberT$resulPaeRasterT))){
          cat(c(paste0(c('grid_', i), collapse = ''), rownames(subset(resulPaeRasterT,
                                                                      resulPaeRasterT[,'grid_n'] == i)), '\n'), sep = '\t', file = logfile, append = TRUE)
        }
      }
    }
    
    # ##################################
    # ## Changing matTemp:
    ######################################################
    matTemp <- preabsMat[, -subtrair] # temporary matrix with homoplastic synapomorphies
    matTemp <- as.matrix(matTemp)
    #######################################################
    
    
    nomesCOL <- colnames(matTemp) # species which are not synapomorphies
    
    print(paste0("Remaining species (non-synapomorphic): ", length(nomesCOL)))
    if(length(nomesCOL) > 0) {
      print(paste0("Species: ", paste(nomesCOL, collapse = ", ")))
    }
    
    quantum <- 1
    for(j in 1:length(nomesCOL)){
      k <- nomesCOL[quantum]
      discre <- matTemp[, as.character(k)]
      for(a in 1:length(discre)){
        ifelse(discre[a] == 1, discre[a] <- k, discre[a] <- 'absent')
      } # vector with the presence and absence of one species
      
      if(all(discre == 'absent')){
        quantum <- quantum + 1
        print(paste0(k, 'is outside of map! Please check it!'))
        break
      }
      
      discre <- as.factor(discre)  ### species presence vector
      
      lista_T[[k]] <- names(discre)[which(discre == k)]
      
      quantum <- quantum + 1
    }
    
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
    
    if(dim(tabela)[1] == 0){
      print('No more iterations are needed.')
      deve_continuar <- FALSE
      break
    }
    
    # Data frame with the results...
    frameTemp <- data.frame(spp = tabela[, 1], grid_n = tabela[, 2], row.names = 1:dim(tabela)[1])
    
    # Matrix with species names and grid numbers:
    resulPaeRaster <- matrix(as.matrix(frameTemp[,2]), nrow(frameTemp),
                             1, dimnames = list(frameTemp[,1], colnames(frameTemp)[2]))
    
    
    ## The number of species supporting the grid number
    speciesNumber <- data.frame(table(resulPaeRaster))
    
    n_occur <- data.frame(table(unlist(lista_T))) # frequency of each grid
    
    
    syn_grids <- unlist(lista_T)[unlist(lista_T) %in% n_occur$Var1[n_occur$Freq > 1]]
    if(length(syn_grids) < 2){
      print('This iteration has to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      next
    }
    
    # ESSENTIAL: Save species lists for each iteration
    lista[[contagem]] <- frameTemp  # species which are not synapomorphies (homoplastic_species)
    listaR[[contagem]] <- frameTempT # species which are synapomorphies (nonHomoplastic_species)
    
    
    preabsMat <- matTemp
    
    print(paste0("New preabsMat matrix: ", dim(preabsMat)[2], " species, ", dim(preabsMat)[1], " grids"))
    
    # FINAL VERIFICATION: If there are not enough species to continue
    if(dim(preabsMat)[2] < 2) {
      print("Less than 2 species remaining - stopping analysis")
      deve_continuar <- FALSE
      break
    }
    
    print(paste0(c('The iteration number ', contagem, ' has finished.'), collapse = ''))
    
    contagem <- contagem + 1
    conta[contagem] <- contagem
    
    print(paste0("=== END OF ITERATION ", contagem - 1, " ==="))
    print("")
  } # close while looping
  
  print("=== WHILE LOOP FINALIZED ===")
  print(paste0("Reason for exit: deve_continuar = ", deve_continuar))
  print(paste0("Characters with CI=1: ", length(which((ciVec == 1) == TRUE))))
  print(paste0("Columns in preabsMat: ", dim(preabsMat)[2]))
  
  # ELEGANT TERMINATIONS FOR END OF LOOP
  if (exists("matTemp") == FALSE && length(conta) == 1){
    # ELEGANT STOP: Temporary matrix does not exist
    return(finalizar_analise_elegante(
      reason = "Temporary matrix was not created",
      details = c("Variable matTemp does not exist after first iteration",
                  "May indicate failure in analysis initialization"),
      iteration = length(conta)
    ))
  }
  
  if(length(colnames(matTemp)) == 0 && length(conta) == 1){
    # ELEGANT STOP: Empty matrix in first iteration
    return(finalizar_analise_elegante(
      reason = "Data matrix empty after first iteration",
      details = c("No columns remaining in matrix after first iteration",
                  "All characters were removed or there is no valid data"),
      iteration = length(conta)
    ))
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
    text(map.r[,c(1, 2)], labels = map.r$gridNumber, cex = 0.8, col = 'black', font = 2)
  }
  
  conta <- NULL
  if(sobrepo == FALSE){
    xis <- seq(1:length(listaR))
    legend(x = 'topright', legend = xis, pch = 15, col = colores[xis],
           title = 'Adding generalized track of the iterations...', title.col = 'red', pt.cex = 1.5, cex = 0.8)
  } else if(sobrepo == TRUE){
    xis <- seq(1:length(listaR))
    legend(x = 'topright', legend = xis, pch = 15, col = colores[xis],
           title = 'Adding generalized track of the iterations...', title.col = 'red', pt.cex = 1.5, cex = 0.8)
  }
  
  if(contagem == 1){
    # ELEGANT TERMINATION: No results
    return(finalizar_analise_elegante(
      reason = "Analysis completed without results",
      details = c("No iteration produced valid results",
                  "Data may not be suitable for PAE-PCE analysis"),
      iteration = contagem,
      plot_grid = FALSE
    ))
  } else if(contagem != 1){
    # NORMAL TERMINATION: Successful analysis
    print("")
    print("")
    print(paste(rep("=", 60), collapse=""))
    print("        PAE-PCE ANALYSIS COMPLETED SUCCESSFULLY")
    print(paste(rep("=", 60), collapse=""))
    print(paste0("Total iterations: ", contagem - 1))
    print(paste(rep("=", 60), collapse=""))
    
    return(list(homoplastic_species = lista,
                nonHomoplastic_species = listaR))
  }
} # close if... else... condition