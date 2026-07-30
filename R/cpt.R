##' Simulated clonal progeny trial
##'
##' This dataset represents a clonal progeny eucalyptus trial. 
##' It contains a simulated phenotype (growth trait, `pheno`) of 600 
##' clones, sampled from 60 families obtained after crossing 20 founders in a 
##' partial diallel design. The trial was laid out 
##' in randomized complete blocks design, with five replicates of single-tree plots. 
##' The coordinates (`row` and `column`) of these plots are also presented. 
##' In the row and column directions, the plants are spaced 
##' 2 m and 3 m apart, respectively.
##' 
##' @docType data
##' 
##' @keywords dataset, simulated
##' 
##' @format ## `cpt`
##'  A data frame with 3000 rows and 6 columns:
##'  \describe{
##'    \item{prog}{60 families}
##'    \item{id}{600 clonal progenies}
##'    \item{row}{40 rows}
##'    \item{col}{75 columns}
##'    \item{block}{5 blocks}
##'    \item{pheno}{3000 phenotypic records}
##'  }
"cpt"
