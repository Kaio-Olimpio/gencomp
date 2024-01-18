##' Eucalyptus real dataset
##'
##' This dataset belongs to clonal trials of the Brazilian pulp company CENIBRA 
##' (Celulose Nipo-Brasileira). It contains the diameter at breast height (`dbh`), 
##' height (`hei`), wood volume (`vol`) and mean annual increment (`mai`) of 100 
##' eucalyptus clones (95 selection candidates and 5 checks). The trial was laid out 
##' in randomized complete blocks design, with 21 replicates and one tree per plot. 
##' This trial do not have contigous blocks: blocks 1-7 are in area 1, blocks 8-14 are in
##' area 2, and blocks 15-21 are in area 3. The `area` columns contains this information.
##' The dataset also has information about the position of each tree in the grid (`row` and `col`). 
##' In the row (dist_row) and column (dist_col) directions, the plants are spaced 
##' 2 m and 3 m apart, respectively. Data were collected at 3 and 6 years after 
##' the trial implementation (information contained in the column `age`).
##' 
##' All rights including Intellectual Property Rights in the Data are owned by CENIBRA and the 
##' authors of this library. Users may only use copies of the Data for the purposes of research, 
##' teaching and education, and for preparing publications, as long as Users acknowledge and 
##' reference the use of this Data in all publications using the citation provided by
##' `citation('competition')`.
##' 
##' @docType data
##' 
##' @keywords dataset
##' 
##' @format ## `eucalyptus`
##'  A data frame with 4200 rows and 13 columns:
##'  \describe{
##'    \item{area}{3 areas}
##'    \item{age}{2 ages}
##'    \item{tree}{2100 trees}
##'    \item{block}{21 blocks}
##'    \item{clone}{100 clones}
##'    \item{dist_row}{2 meters between rows}
##'    \item{dist_col}{3 meters between columns}
##'    \item{row}{42 rows}
##'    \item{col}{50 columns}
##'    \item{dbh}{Diameter at breast heigh (cm)}
##'    \item{hei}{Height (m)}
##'    \item{vol}{Wood volume (m<sup>3</sup>)}
##'    \item{mai}{Mean annual increment of wood volume (m<sup>3</sup> ha<sup>-1</sup> year<sup>-1</sup>)}
##'  }
"eucalyptus"
