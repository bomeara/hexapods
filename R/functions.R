gettrees <- function(clade="Hexapoda") {
  treezip <- "http://files.opentreeoflife.org/synthesis/opentree9.1/opentree9.1_tree.tgz"
  download.file(treezip,destfile="opentree9.1_tree.tgz")
  untar("opentree9.1_tree.tgz")
  # supertree is taxonomy plus source trees; grafted is only source trees
  supertree.phy <- ape::read.tree("opentree9.1_tree/labelled_supertree/labelled_supertree_ottnames.tre")
  grafted.phy <- ape::read.tree("opentree9.1_tree/grafted_solution/grafted_solution_ottnames.tre")
  supertree.phy <- ape::extract.clade(supertree.phy, supertree.phy$node.label[which(grepl(clade, supertree.phy$node.label))])
  grafted.phy <- ape::extract.clade(grafted.phy, grafted.phy$node.label[which(grepl(clade, grafted.phy$node.label))])
  return(list(supertree.phy=supertree.phy, grafted.phy=grafted.phy ))
}
