# plan <- drake_plan(
#   trees = get_trees(),
#   ot_taxonomy = get_ot_taxonomy(),
#   print(trees),
#   save(trees, file=file_out("ot_rawtrees.rda")),
#   ot_supertree = convert_tiplabels_to_genbank_fast_parse(trees$supertree.phy),
#   save(ot_supertree, file=file_out("ot_supertree.rda")),
#   ot_graftedtree = convert_tiplabels_to_genbank_fast_parse(trees$grafted.phy),
#   save(ot_graftedtree, file=file_out("ot_graftedtree.rda")),
#   taxa_from_tree = extract_names(ot_supertree),
#   save(taxa_from_tree, file=file_out("taxa_from_tree.rda")),
#   #taxa_with_genbank = wrap_seqs_in_genbank(taxa_from_tree),
#   #save(taxa_with_genbank, file=file_out("taxa_with_genbank.rda")),
#   #taxa_with_funding = get_funding(taxa_with_genbank),
#   taxa_with_funding = get_funding(taxa_from_tree),
#   save(taxa_with_funding, file=file_out("taxa_with_funding.rda")),
#   report = rmarkdown::render(
#     knitr_in("report.Rmd"),
#     output_file = file_out("report.html"),
#     quiet = TRUE
#   )
# )

# familyrun <- drake_plan(
#   families = get_families(),
#   taxa_with_funding = get_funding(data.frame(taxon=families, stringsAsFactors=FALSE)),
#   write.csv(taxa_with_funding, file=file_out("taxa_with_funding.csv")),
#   taxa_with_funding_genbank = wrap_seqs_in_genbank(taxa_with_funding),
#   write.csv(taxa_with_funding_genbank, file=file_out("taxa_with_funding_genbank.csv")),
#   taxa_with_funding_genbank_dark = wrap_dark_in_genbank(taxa_with_funding_genbank),
#   write.csv(taxa_with_funding_genbank_dark, file=file_out("taxa_with_funding_genbank_dark.csv"))
# )


familyrun <- drake_plan(
  taxa_df = get_hexapoda_info(),
  write.csv(taxa_df, file=file_out("docs/raw_taxa.csv"))
)
#   taxa_genbank = wrap_seqs_in_genbank(taxa_df),
#   write.csv(taxa_with_funding_genbank, file=file_out("taxa_with_funding_genbank.csv")),
#   taxa_with_funding_genbank_dark = wrap_dark_in_genbank(taxa_with_funding_genbank),
#   write.csv(taxa_with_funding_genbank_dark, file=file_out("taxa_with_funding_genbank_dark.csv"))
# )


taxarun <- drake_plan(
  taxon.names = extract_taxon_info_from_dir_of_papers(),
  write.csv(taxon.names, file=file_out("taxa_from_papers.csv"))
)
#
# scholarplan <- drake_plan(
#   scholarinfo <- loop_counts_from_scholar("Coleoptera_fams_subfams_01Aug2018_NPL")
#   write.csv(scholarinfo, file="ScholarCounts.csv")
# )
