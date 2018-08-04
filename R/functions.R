# see also https://raw.githubusercontent.com/bomeara/phyloview/master/phyloview.Rmd and rphylotastic's darktaxa section

get_genbank_taxonomy <- function(clade="Hexapoda") {
  # TO DO : download genbank taxonomy file, so can match up ncbi ids from ot with genbank names
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

wrap_dark_in_genbank <- function(taxon.dataframe) {
  taxon.dataframe$genbank.dark.count <- NA
  taxon.dataframe$genbank.known.count <- NA
  taxon.dataframe$genbank.dark.fraction <- NA
  for (i in sequence(nrow(taxon.dataframe))) {
    local.df <- dark_in_genbank(taxon.dataframe$taxon[i])
    if(nrow(local.df)==1) {
      taxon.dataframe$genbank.dark.count[i] <- local.df$dark.count
      taxon.dataframe$genbank.known.count[i] <- local.df$known.count
      taxon.dataframe$genbank.dark.fraction[i] <- local.df$fraction.dark
    }
  }
  return(taxon.dataframe)
}

dark_in_genbank <- function(taxon) {
  result <- rphylotastic::taxon_separate_dark_taxa_using_genbank(taxon, sleep=3)
  result.df <- data.frame(dark.count=length(result$dark), known.count=length(result$known), fraction.dark=result$fraction.dark)
  return(result.df)
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

get_funding <- function(taxon.dataframe, stem_name=TRUE) {
  data(grants) #from rnsf
  good.grant.indices <- which(grepl("systematics|phylo|bioinfor|taxonom|revision|peet", grants$fundProgramName, ignore.case=TRUE))
  good.grant.indices <- c(good.grant.indices,which(grepl("systematics|phylo|bioinfor|taxonom|revision|peet", grants$abstractText, ignore.case=TRUE)))
  good.grant.indices <- c(good.grant.indices,which(grepl("systematics|phylo|bioinfor|taxonom|revision|peet", grants$title, ignore.case=TRUE)))


  relevant.grants <- grants[unique(good.grant.indices),]
  get_funding_for_taxon <- function(taxon, r.grants, stem_name=TRUE) {
    funding <- 0
    if(stem_name) {
      taxon <- gsub("dae", "d", taxon)
    }
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
  #TODO Add '+OR+' between multiple entries (which will come in as family)
  page <- xml2::read_html(paste0('https://scholar.google.com/scholar?hl=en&as_sdt=0%2C14&q=allintitle%3A+', family, '+topology+OR+phylogeny+OR+phylogenetics+OR+cladistics+OR+phylogenetic+OR+cladistic&btnG='))
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
      Sys.sleep(60*15)
    }
  }
  return(taxonsheet)
}

get_figures_from_plosone <- function(query="phylogeny", dir=NULL) {
  original.dir <- getwd()
  if(!is.null(dir)) {
    if(!dir.exists(dir)) {
      dir.create(dir)
    }
    setwd(dir)
  }
  results <- rplos::searchplos(q= query, fl= "id", limit = 1000, fq='journal_key:PLoSONE')$data$id
  for (i in seq_along(results)) {
    for (fig.index in sequence(2)) {
      try(utils::download.file(paste0('http://journals.plos.org/plosone/article/figure/image?size=large&download=&id=', results[i], '.g00', fig.index), destfile=paste0('Fig_', gsub('/','_',results[i]), '.g00', fig.index, '.png'), method="internal"))
      Sys.sleep(6)
    }
  }
  setwd(original.dir)
}

train_tree_model <- function(train_dir = "training", validation_dir="validation") {

  # model <- application_vgg16(
  #   weights = "imagenet",
  #   include_top = FALSE
  # )

  test_datagen <- image_data_generator(rescale = 1/255)

  train_generator <- flow_images_from_directory(
  train_dir,                  # Target directory
  test_datagen,              # Data generator
  target_size = c(150, 150),  # Resizes all images to 150 × 150
  batch_size = 20,
  class_mode = "binary"       # binary_crossentropy loss for binary labels
  )

  validation_generator <- flow_images_from_directory(
    validation_dir,
    test_datagen,
    target_size = c(150, 150),  # Resizes all images to 150 × 150
    batch_size = 20,
    class_mode = "binary"
  )

#   model %>% compile(
#   loss = "binary_crossentropy",
#   optimizer = optimizer_rmsprop(lr = 2e-5),
#   metrics = c("accuracy")
# )

  conv_base <- application_vgg16(
    weights = "imagenet",
    include_top = FALSE,
    input_shape = c(150, 150, 3)
  )

  model <- keras_model_sequential() %>%
    conv_base %>%
    layer_flatten() %>%
    layer_dense(units = 256, activation = "relu") %>%
    layer_dense(units = 1, activation = "sigmoid")

    freeze_weights(conv_base)

  opt <- optimizer_rmsprop(lr = 0.0001, decay = 1e-6)

  model %>% compile(
    loss = "mse",
    optimizer = opt,
    metrics = "accuracy"
  )

  history <- model %>% fit_generator(
   train_generator,
   steps_per_epoch = 30,
   epochs = 50,
   validation_data = validation_generator,
   validation_steps = 3
  )
  return(model)
}

predict_using_model <- function(dir, model) {
  original.dir <- getwd()
  setwd(dir)
  image_prep <- function(x) {
    arrays <- lapply(x, function(dir) {
      img <- image_load(dir, target_size = c(150,150))
      x <- image_to_array(img)
      x <- array_reshape(x, c(1, dim(x)))
      x <- imagenet_preprocess_input(x)
    })
    do.call(abind::abind, c(arrays, list(along = 1)))
  }
  files <- system(paste0("ls -1 ", dir), intern=TRUE)
  res <- predict(model, image_prep(files))[,1]
  names(res) <- files
  setwd(original.dir)
  return(res)
}

#' Runs on one line
extract_classification_from_gnr_resolve <- function(x) {
  result <- data.frame(t(strsplit(x['classification_path'], "\\|",)[[1]]), stringsAsFactors=FALSE)
  colnames(result) <- t(strsplit(x['classification_path_ranks'], "\\|",)[[1]])
  result <- cbind(data.frame(taxon=x['matched_name'], user_supplied_name=x['user_supplied_name'], taxon_rank=colnames(result)[ncol(result)], stringsAsFactors=FALSE), result)
  return(result)
}

extract_taxon_info_from_paper <- function(file) {
  taxa <- rphylotastic::file_get_scientific_names(file)
  taxa.resolved <-  as.data.frame(taxize::gnr_resolve(taxa, data_source_ids = 1, best_match_only = TRUE, fields="all"), stringsAsFactors=FALSE)
  all.taxa <- plyr::rbind.fill(apply(taxa.resolved, 1,extract_classification_from_gnr_resolve))
}

extract_taxon_info_from_dir_of_papers <- function(path='/Users/bomeara/Google Drive/SoToL Hexapods (File responses)/Full paper (File responses)/') {
  original.dir <- getwd()
  setwd(path)
  files <- system("ls -1", intern=TRUE)
  results <- data.frame()
  for (i in seq_along(files)) {
    local.results <- data.frame()
    try(local.results <- extract_taxon_info_from_paper(files[i]))
    if(nrow(local.results)>0) {
      local.results$paper <- files[i]
      results <- plyr::rbind.fill(results, local.results)
    }
  }
  setwd(original.dir)
  return(results)
}

get_families <- function() {
  # insecta: ae304a1e0beadcfec04932589049bb5a
  families_raw <- taxize::downstream('ae304a1e0beadcfec04932589049bb5a', downto = 'Family', db = 'col')[[1]]$childtaxa_name
  return(families_raw[!grepl(' ', families_raw)] ) #get rid of ones that are "Not assigned"
}

get_otol_tree <- function(taxa) {
  taxa <- as.character(taxa)
  ott_ids <- rep(NA, length(taxa))
  taxa.chunks <- split(taxa, ceiling(seq_along(taxa)/250))
  position <- 0
  for (i in sequence(length(taxa.chunks))) {
    ott_ids[(position+1):(position+length(taxa.chunks[[i]]))] <- rotl::tnrs_match_names(taxa.chunks[[i]])$ott_id
    position <- position + length(taxa.chunks[[i]])
  }
  ott_ids <- ott_ids[!is.na(ott_ids)]

}
