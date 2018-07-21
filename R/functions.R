# see also https://raw.githubusercontent.com/bomeara/phyloview/master/phyloview.Rmd and rphylotastic's darktaxa section

get_trees <- function(clade="Hexapoda") {
  treezip <- "http://files.opentreeoflife.org/synthesis/opentree9.1/opentree9.1_tree.tgz"
  download.file(treezip,destfile="opentree9.1_tree.tgz")
  untar("opentree9.1_tree.tgz")
  # supertree is taxonomy plus source trees; grafted is only source trees
  supertree.phy <- ape::read.tree("opentree9.1_tree/labelled_supertree/labelled_supertree_ottnames.tre")
  grafted.phy <- ape::read.tree("opentree9.1_tree/grafted_solution/grafted_solution_ottnames.tre")
  supertree.phy <- ape::extract.clade(supertree.phy, supertree.phy$node.label[which(grepl(clade, supertree.phy$node.label))])
  grafted.phy <- ape::extract.clade(grafted.phy, grafted.phy$node.label[which(grepl(clade, grafted.phy$node.label))])

  #for debugging
  grafted.phy <- geiger::drop.random(grafted.phy, ape::Ntip(grafted.phy)-25)
  supertree.phy <- geiger::drop.random(supertree.phy, ape::Ntip(supertree.phy)-50)


  return(list(supertree.phy=supertree.phy, grafted.phy=grafted.phy ))
}

convert_tiplabels_to_genbank <- function(phy) {
  try_names <- function(x) {
    if(nchar(x)==0) {
      return(x)
    }
    final.name <- x
    ott.id <- gsub(".*_ott","",x)
    new.final <- NA
    try(new.final <- rotl::tax_name(rotl::taxonomy_taxon_info(ott.id))[[1]])
    if(!is.na(new.final)) {
      final.name <- new.final
    }
    return(final.name)
  }
  phy$tip.label <- sapply(phy$tip.label, try_names)
  phy$node.label <- sapply(phy$node.label, try_names)
  return(phy)
}

wrap_seqs_in_genbank <- function(taxon.dataframe) {
  taxon.dataframe$genbank.count <- sapply(taxon.dataframe$taxon, seqs_in_genbank)
  return(taxon.dataframe)
}

seqs_in_genbank <- function(taxon, genes=c("COI", "cytochrome", "18S", "28S", "enolase", "elongation factor"), sleep.time=3) {
  #sitophilus oryzae[organism] AND (COI Or Cytochrome OR 28S OR 18S OR enolase OR elongation)
  query <- paste0(taxon, '[organism] AND (',paste(genes, collapse=" OR "),')')
  Sys.sleep(sleep.time) # to keep from clobbering NCBI
  result <- NA
  try(result <- rentrez::entrez_search(db="nuccore", query, use_history=TRUE)$count)
  return(result)
}

extract_names <- function(phy) {
  all.names <- c(phy$tip.label, phy$node.label)
  node.numbers <- sequence(length(all.names))
  missing.names <- which(nchar(all.names)==0)
  all.names <- all.names[-missing.names]
  node.numbers <- node.numbers[-missing.names]
  final.df <- data.frame(taxon=all.names, node.id = node.numbers, is.tip = FALSE, stringsAsFactors=FALSE)
  final.df$is.tip[final.df$node.id<=ape::Ntip(phy)] <- TRUE
  return(final.df)
}
