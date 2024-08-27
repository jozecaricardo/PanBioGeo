### Function to eliminate singletons and return a data frame object ###
singleton.to.data.frame <- function(data = NULL, phylogeny = NULL){
	lat <- data[,2]
	long <- data[,3]
	names(lat) <- data[, 1]
	names(long) <- data[, 1]
	locais <- cbind(lat, long)

	# taxon.tree <- read.tree(file = phylogeny, collapse = NULL)
	taxon.tree <- phylogeny
	
	singletons <- NULL # apenas uma ocorrência
	sppsORD <- sort(rownames(locais), decreasing = F)
	for(i in 1:length(unique(rownames(locais)))){
  		singletons[i] <- length(grep(pattern = unique(sppsORD)[i], x = sppsORD))
	}

	singletons <- which(singletons == 1)
	um_ponto <- unique(sppsORD)[singletons]

	# Modifying the species tree:
	taxonNew.tree <- drop.tip(taxon.tree, um_ponto)
	# plotTree(plebeNew.tree)

	retira <- name.check(taxonNew.tree, locais) #olhando os nomes na arvore e nos dados...
	
	# now without singletons:
	locaisNew <- locais
	for(i in 1:length(unique(retira$data_not_tree))){
  		linha <- which(rownames(locaisNew) == unique(retira$data_not_tree)[i])
  		locaisNew <- locaisNew[-linha,]
	}

	# data frame
	lat <- locaisNew[, 1]
	long <- locaisNew[, 2]
	spp <- rownames(locaisNew)

	data_df <- data.frame(spp = spp, long = long, lat = lat)

	return(list(data_df = data_df, treeMod = taxonNew.tree))
}