##' Eucalyptus dataset
##'
##' This dataset represents a clonal eucalyptus trial. It contains the mean annual increment (`MAI`) of 100 
##' clones. The trial was laid out 
##' in randomized complete blocks design, with 13 replicates and one tree per plot. 
##' This trial did not have contiguous blocks: the first six blocks were in area 1 and the last seven were in
##' area 2. The `area` column contains this information.
##' The dataset also has information about the position of each tree in the grid (`row` and `col`). 
##' In the row and column directions, the plants are spaced 
##' 2 m (dist_row) and 3 m (dist_col) apart, respectively. Data were collected at 3 and 6 years after 
##' the trial implementation (information contained in the column `age`).
##' 
##' 
##' @docType data
##' 
##' @keywords dataset
##' 
##' @format ## `euca`
##'  A data frame with 4200 rows and 10 columns:
##'  \describe{
##'    \item{area}{2 areas}
##'    \item{age}{2 ages}
##'    \item{tree}{1144 trees}
##'    \item{block}{13 blocks}
##'    \item{clone}{100 clones}
##'    \item{dist_row}{2 meters between rows}
##'    \item{dist_col}{3 meters between columns}
##'    \item{row}{26 rows}
##'    \item{col}{44 columns}
##'    \item{MAI}{Mean annual increment of wood volume (m<sup>3</sup> ha<sup>-1</sup> year<sup>-1</sup>)}
##'  }
"euca"
