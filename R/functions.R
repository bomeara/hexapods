# see also https://raw.githubusercontent.com/bomeara/phyloview/master/phyloview.Rmd and rphylotastic's darktaxa section

get_genbank_taxonomy <- function(clade="Hexapoda") {
  # TO DO : download genbank taxonmy file, so can match up ncbi ids from ot with genbank names
}

get_ot_taxonomy <- function(clade="Hexapoda") {
  taxzip <- "http://files.opentreeoflife.org/ott/ott3.0/ott3.0.tgz"
  download.file(taxzip,destfile="ott3.0.tgz")
  untar("ott3.0.tgz")
  ott.taxonomy <- readr::read_tsv("ott/taxonomy.tsv")
  ott.taxonomy <- as.data.frame(ott.taxonomy)
  ott.taxonomy <- ott.taxonomy[,!grepl("\\|", colnames(ott.taxonomy))] #dealing with pipes
  ott.taxonomy <- ott.taxonomy[,!grepl("X", colnames(ott.taxonomy))]
  extract_ncbi <- function(x) {
    result <- gsub("ncbi:","",regmatches(x,regexpr("ncbi:\\d+", x, perl=TRUE)))
    if(length(result)==0) {
      result <- NA
    }
    return(result)
  }
  ott.taxonomy$ncbi.number <- sapply(ott.taxonomy$sourceinfo, extract_ncbi)
  return(ott.taxonomy)
}

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
  #grafted.phy <- geiger::drop.random(grafted.phy, ape::Ntip(grafted.phy)-25)
  #supertree.phy <- geiger::drop.random(supertree.phy, ape::Ntip(supertree.phy)-50)

  #grafted.phy <- convert_tiplabels_to_genbank(grafted.phy)
  #supertree.phy <- convert_tiplabels_to_genbank(supertree.phy)


  return(list(supertree.phy=supertree.phy, grafted.phy=grafted.phy ))
}

convert_tiplabels_to_genbank_full_parse <- function(phy) {
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

convert_tiplabels_to_genbank_fast_parse <- function(phy) {
  try_names_fast <- function(x) {
    if(nchar(x)==0) {
      return(x)
    }
    final.name <- strsplit(x,"_ott",x)[[1]][1]
    final.name <- gsub("_", " ", final.name)
    if(nchar(final.name) == nchar(x)) { #taxon Campsicnemussp.2KRG-2014ott5829692 is an example
      ott.id <- gsub(".*ott","",x)
      new.final <- NA
      try(new.final <- rotl::tax_name(rotl::taxonomy_taxon_info(ott.id))[[1]])
      if(!is.na(new.final)) {
        final.name <- new.final
      }
    }
    return(final.name)
  }
  phy$tip.label <- sapply(phy$tip.label, try_names_fast)
  phy$node.label <- sapply(phy$node.label, try_names_fast)
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
  bad.names <- which(nchar(all.names)==0)
  bad.names <- unique(c(which(grepl("mrcaott", all.names))))
  if(length(bad.names)>0) {
    all.names <- all.names[-bad.names]
    node.numbers <- node.numbers[-bad.names]
  }
  final.df <- data.frame(taxon=all.names, node.id = node.numbers, is.tip = FALSE, stringsAsFactors=FALSE)
  final.df$is.tip[final.df$node.id<=ape::Ntip(phy)] <- TRUE
  return(final.df)
}

get_funding <- function(taxon.dataframe) {
  data(grants) #from rnsf
  good.grant.indices <- which(grepl("systematics|phylo|bioinfor|taxonom|revision", grants$fundProgramName, ignore.case=TRUE))
  good.grant.indices <- c(good.grant.indices,which(grepl("systematics|phylo|bioinfor|taxonom|revision", grants$abstractText, ignore.case=TRUE)))
  good.grant.indices <- c(good.grant.indices,which(grepl("systematics|phylo|bioinfor|taxonom|revision", grants$title, ignore.case=TRUE)))


  relevant.grants <- grants[unique(good.grant.indices),]
  get_funding_for_taxon <- function(taxon, r.grants) {
    funding <- 0
    matching.grant.indices <- unique(c(which(grepl(taxon, r.grants$abstractText, ignore.case=TRUE)), which(grepl(taxon, r.grants$title, ignore.case=TRUE))))
    if(length(matching.grant.indices)>0) {
      funding <- sum(as.numeric(r.grants$fundsObligatedAmt[matching.grant.indices]), na.rm=TRUE)
    }
    return(funding)
  }
  taxon.dataframe$nsf.funding <- sapply(taxon.dataframe$taxon, get_funding_for_taxon, r.grants=relevant.grants)
  return(taxon.dataframe)
}

get_counts_from_scholar <- function(family) {
  page <- xml2::read_html(paste0('https://scholar.google.com/scholar?hl=en&as_sdt=0%2C14&q=allintitle%3A+', family, '+OR+',gsub("dae", "d", family),'+topology+OR+phylogeny+OR+phylogenetics+OR+cladistics+OR+phylogenetic+OR+cladistic&btnG='))
  section <- rvest::html_nodes(page, '.gs_ab_mdw')[2]
  result <- gsub(",", "", gsub("About ", "", stringr::str_extract(as.character(section), "About \\d+\\,?\\d*+")))
  return(result)
}

loop_counts_from_scholar <- function(sheet_name) {
  taxonsheet <- googlesheets::gs_read(googlesheets::gs_title(sheet_name))
  taxonsheet$ScholarCount <- 0
  for(i in sequence(nrow(taxonsheet))) {
    if(!is.na(taxonsheet$Family[i])) {
      try(taxonsheet$ScholarCount[i] <- get_counts_from_scholar(taxonsheet$Family[i]))
      print(paste(taxonsheet$Family[i], taxonsheet$ScholarCount[i]))
      Sys.sleep(10)
    }
  }
  return(taxonsheet)
}
