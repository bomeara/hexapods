plan <- drake_plan(
  trees = get_trees(),
  print(trees),
  ot_supertree = convert_tiplabels_to_genbank(trees$supertree.phy),
  ot_graftedtree = convert_tiplabels_to_genbank(trees$grafted.phy),
  save(ot_supertree, file="ot_supertree.rda"),
  taxa_from_tree = extract_names(ot_supertree),
  taxa_with_genbank = wrap_seqs_in_genbank(taxa_from_tree),
  save(taxa_with_genbank, file="taxa_with_genbank.rda"),
  report = rmarkdown::render(
    knitr_in("report.Rmd"),
    output_file = file_out("report.html"),
    quiet = TRUE
  )
)
