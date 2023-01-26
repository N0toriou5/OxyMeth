#### Function to convert a human regulon to any species
require("babelgene")

#### Function to translate a human regulon to another specie
Regulon_HumantoModelSym <- function(regulon, species = "mouse", limit = NULL) {
  
  ### limit to nodes showing targets >= limit
  if (is.null(limit)) {
    
    regulon <- regulon
    
    } else {
      limit <- limit 
      ntar <- sapply(regulon, function(x){
        length(x$tfmode)
      })
      keep <- which(ntar >= limit)
      regulon <- regulon[keep]
    }
  
  gene_mapping <- orthologs(genes = names(regulon), species = species)
  human_genes <- c(names(regulon))
  member_names <- unique(intersect(names(regulon), gene_mapping$human_symbol))
  index <- which(human_genes %in% member_names)
  regulon <- regulon[index]
  idx <- match(names(regulon), gene_mapping$human_symbol)
  gene_mapping <- gene_mapping[idx,]
  #gene_mapping <- gene_mapping[!duplicated(gene_mapping$symbol),]
  member_names <- unique(intersect(names(regulon), gene_mapping$human_symbol))
  index <- which(human_genes %in% member_names)
  
  #### convert centroids
  for (i in 1:nrow(gene_mapping)) {
    human_gene <- gene_mapping[i,1]
    rat_map <- gene_mapping[gene_mapping$human_symbol == human_gene,]
    if (nrow(rat_map) >= 2) {
      idx <- max(rat_map$support_n)
      rat_gene <- rat_map[rat_map$support_n==idx, 5]
      names(regulon)[i] <- rat_gene
    } else {
      rat_gene <- rat_map[rat_map$human_symbol == human_gene, 5]
      names(regulon)[i] <- rat_gene
    }
  }
  
  #### now for the target
  for (i in 1:length(regulon)) {
    nameVec <- names(regulon[[i]]$tfmode)
    tryCatch(gene_mapping <- orthologs(genes = nameVec, species = species),
             error = function(x){         
               message("ERROR: some nodes in the regulon may have 0 targets left.")
               stop("Please, specify a size limit for the sub-networks!")
             })
    idx <- match(nameVec, gene_mapping$human_symbol)
    gene_mapping <- gene_mapping[idx,]
    names(regulon[[i]]$tfmode) <- gene_mapping$symbol
    na.omit(names(regulon[[i]]$tfmode))
    id <- which(is.na(names(regulon[[i]]$tfmode)))
    regulon[[i]]$tfmode <- regulon[[i]]$tfmode[-id]
    regulon[[i]]$likelihood <- regulon[[i]]$likelihood[-id]        
  }
  class(regulon)<-"regulon"
  return(regulon)
}
