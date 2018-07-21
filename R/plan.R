plan <- drake_plan(
  trees = get_trees(),
  print(trees),
  ot_supertree = convert_tiplabels_to_genbank(trees$supertree.phy),
  ot_graftedtree = convert_tiplabels_to_genbank(trees$grafted.phy,
  report = rmarkdown::render(
    knitr_in("report.Rmd"),
    output_file = file_out("report.html"),
    quiet = TRUE
  )
)
