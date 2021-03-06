---
title: "Analyses from the FuturePhy Hexapod group"
author: Brian O'Meara
output: html_document
---

# State of the hexapod tree of life

This is work *in progress* from the Hexapod FuturePhy group (meeting in Chicago, Aug. 1-3, 2018). It will update as we add information.



For example, we are planning to gather info on 458 families so far, but other families will be added.

The code for these analyses is all in https://github.com/bomeara/hexapods. The ROpenSci package `drake` is used to dynamically run analyses to generate content (including this page).


## Binomials per clade in GenBank

![plot of chunk plotgenbankknown](figure/plotgenbankknown-1.png)

## Dark taxa per clade in GenBank

![plot of chunk plotgenbankdark](figure/plotgenbankdark-1.png)

## BOLD taxa per clade in GenBank

![plot of chunk plotgenbankbold](figure/plotgenbankbold-1.png)

## Species in GenBank vs diversity in Catalog of Life

![plot of chunk plotgbvscol](figure/plotgbvscol-1.png)

Line shows 100% of Catalog of Life species in GenBank (note that this is based on species counts in families, not 1:1 matching). File with info on each family is [taxa_col.csv.gz](taxa_col.csv.gz); file with info on each species from the Catalog of Life is [all_species_hierarchy.csv.gz](all_species_hierarchy.csv.gz). See [here](https://github.com/bomeara/hexapods/blob/master/R/plan.R) for the R code to organize data generation and [here](https://github.com/bomeara/hexapods/blob/master/R/functions.R) for the functions used. Packages used are listed [here](https://github.com/bomeara/hexapods/blob/master/R/packages.R): main ones are `taxize`, `rphylotastic`, `googlesheets`, `rotl`, and `rentrez`.

## Funding

The core grant for this was NSF [FuturePhy grant ID]. Participation and tools was also supported by grants [GN], [Phylotastic], [more]
