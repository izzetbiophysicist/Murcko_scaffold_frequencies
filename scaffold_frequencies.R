library(rcdk)

args = commandArgs(trailingOnly=TRUE)

all_compounds <- read.csv(args[1], row.names = 1)

count_scaffolds <-function(smiles_list){
  sclist <- c()
  for(i in 1:length(smiles_list)){
  
    m1 <- parse.smiles(smiles_list[i])
    a <-  tryCatch(get.murcko.fragments(m1, single.framework = TRUE), error=function(e) NULL)
    sclist <- c(sclist, a[[1]]$frameworks)
  
  }
  scunique <- unique(sort(sclist))
  
  if(length(scunique) > 0){
  scaf_data <- data.frame(scaffold=scunique, count=0)
  for(q in 1:nrow(scaf_data)){
    scaf_data[q,2] <- length(which(sclist == scaf_data[q,1]))
  }
  return(scaf_data)
  }else{
    return(data.frame(scaffold=NA, count=NA))
  }
}

scaffolds <- count_scaffolds(all_compounds$smiles)

write.csv(scaffolds, paste(args[2], "scaffolds.csv", sep=""))