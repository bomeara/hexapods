plan <- drake_plan(
  trees = get_trees(),
  print(trees),
  save(trees, file=file_out("ot_rawtrees.rda")),
  ot_supertree = convert_tiplabels_to_genbank_fast_parse(trees$supertree.phy),
  save(ot_supertree, file=file_out("ot_supertree.rda")),
  ot_graftedtree = convert_tiplabels_to_genbank_fast_parse(trees$grafted.phy),
  save(ot_graftedtree, file=file_out("ot_graftedtree.rda")),
  taxa_from_tree = extract_names(ot_supertree),
  save(taxa_from_tree, file=file_out("taxa_from_tree.rda")),
  taxa_with_genbank = wrap_seqs_in_genbank(taxa_from_tree),
  save(taxa_with_genbank, file=file_out("taxa_with_genbank.rda")),
  taxa_with_gb_and_funding = get_funding(taxa_with_genbank),
  save(taxa_with_gb_and_funding, file=file_out("taxa_with_gb_and_funding.rda")),
  report = rmarkdown::render(
    knitr_in("report.Rmd"),
    output_file = file_out("report.html"),
    quiet = TRUE
  )
)
