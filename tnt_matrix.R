tnt_matrix <- function(pres_abs){
#pres_abs = .txt file or an object which is the output from the node_terminal function
  if(class(pres_abs)[1] == 'matrix' || class(pres_abs)[2] == 'array'){
    mat <- pres_abs
  } else {
    mat <- read.table(file = pres_abs, sep = '', dec = '.', row.names = 1) # presence-absence matrix
  }

  #extraindo os número das células do raster
  taxa <- row.names(mat)

  #produzindo os arquivos *.txt
  n.cara <- ncol(mat)
  n.taxa <- nrow(mat)

  nomes <- NULL

  conta <- 0

  # capture R output:
  zz <- textConnection('texto', "w")
  sink(zz)
  cat('nstates num 32;',sep = '\n')
  cat('xread',sep = '\n')
  #cat('/* Matrix for a PAE-PCE analysis */', sep = '\n')
  cat(n.cara, n.taxa, fill = T)

  # ROOT:
  # cat(paste0(c('ROOT', ' ', rep(0, n.cara)), collapse = ''), sep = '\n')

  #colocando a matriz de caracteres pres_abs...
  for(j in 1:n.taxa){
    nomes[j] <- paste(c(taxa[j]),collapse = ', ')
    cat(paste(mat[j,], collapse = ''), fill = T,
      labels = paste0(nomes[j]))
  }

  cat(';',sep = '\n')

  # taxon names:
  cat(sep = '\n')
  cat('cnames', sep = '\n')
  for(i in 1:n.cara){
    cat(paste0(c('{',conta, ' ', colnames(mat)[i], ';'), collapse = ''), sep = '\n')
    conta <- conta + 1
  }

  cat(';',sep = '\n')
  cat(sep = '\n')
  # cat('cc+ 0.6;',sep = '\n')
  cat('proc/;',sep = '\n')
  sink()
  close(zz)

  #abrindo, escrevendo e fechando o arquivo...
  exte <- 'out/pres_abs.tnt'
  tempor = file(exte, "wt")
  writeLines(texto, con = tempor)
  close(tempor)
}