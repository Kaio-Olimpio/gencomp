##' @title Extract outputs from the genetic-spatial competition model
##' 
##' @description
##' This function provides responses extracted from the genetic-spatial 
##' competition model fitted using [gencomp::asr()]. 
##' 
##' @param prep.out A `comprep` object.
##' @param model An `asreml` object, preferably obtained using the [gencomp::asr()] function.
##' @param weight.tgv A logical value. If `TRUE`, the function will use the direct and
##' indirect genetic effects' reliability as a weight when estimating the total genotypic
##' value. Defaults to `FALSE`. See Details for more information on these methods.
##'
##' @return The function returns:
##' \itemize{
##' \item \code{lrt} : a data frame with the Likelihood Ratio Test results.
##' \item \code{varcomp} : A data frame summarizing the random parameter vector 
##' (object$vparameters). Variance component ratios are included if param = "gamma", 
##' and a measure of precision (standard error) is included along with boundary 
##' constraints at termination and the percentage change in the final iteration.
##' \item \code{blups} : A list with a dataframe containing the direct (DGE) and indirect genetic effects (IGE), their standard errors, 
##' the competition class of each genotype and the total genotypic value (TGV). If `age = TRUE` in 
##' [gencomp::prep()], `blup` will contain a further dataframe the within-ages DGE and IGE. If 
##' other random effects were declared in the model, `blup` will contain a last data frame 
##' with their BLUPs.
##' }
##'
##' @details
##' The genetic-spatial competition model provides the direct (DGE) and indirect 
##' genetic effects (IGE) of each genotype. The DGE represents the "pure" performance
##' of the genotype, while the IGE is the related to the average effect of the genotype on
##' the genotypic value of its neighbours. The higher the IGE, the more aggressive is the genotype. 
##' Here, we use the classification proposed by \insertCite{ferreira_novel_2023;textual}{gencomp} 
##' to define competition classes: 
##' 
##' \deqn{\begin{cases} c_i > \overline{c} + sd(c) \rightarrow \text{Aggressive} \\ \overline{c} + sd(c) > c_i > \overline{c} - sd(c) \rightarrow \text{Homeostatic} \\ c_i < \overline{c} - sd(c) \rightarrow \text{Sensitive} \end{cases}}
##' 
##' where \eqn{c_i} is the IGE of the i<sup>th</sup> genotype, \eqn{\overline{c}} is 
##' the mean IGE, and \eqn{sd(c)} is the standard deviation of the IGE. This classification is 
##' detailed in the plot `IGE.density`.
##' 
##' The total genotypic value (TGV) is given by: 
##' 
##' \deqn{TGV_i =  d_i + CIF \times c_i}
##' 
##' where \eqn{d_i} is the DGE of the i<sup>th</sup> genotype, and \eqn{CIF} is 
##' the overall competition intensity factor, previously computed in the function 
##' [gencomp::prep()]. If `weight.tgv = TRUE`, the DGE and IGE will be 
##' multiplied by their respective reliabilities (\eqn{r_{d_i}^2} and \eqn{r_{c_i}^2}):
##' 
##' \deqn{wTGV_i = r_{d_i}^2 \times d_i + r_{c_i}^2 \times {CIF \times c_i}}
##'
##' @seealso  [gencomp::prep], [gencomp::asr]
##' 
##' @references 
##' \insertAllCited{}
##' 
##' @import ggplot2
##' @importFrom rlang .data
##' @importFrom stats reorder sd density
##' @importFrom Rdpack reprompt
##' @export
##' 
##' @examples
##' \donttest{
##'  library(gencomp)
##'  comp_mat = prep(data = euca, gen = 'clone', repl = 'block', area = 'area', 
##'                  ind = 'tree', age = 'age', row = 'row', col = 'col', 
##'                  dist.col = 3, dist.row = 2, trait = 'mai', method = 'SK',
##'                  n.dec = 3, verbose = TRUE)
##'  model = asr(prep.out = comp_mat, 
##'              fixed = mai~ age, 
##'              random = ~ block:age, 
##'              cor = TRUE, maxit = 20,
##'              lrtest = FALSE)
##'                   
##'  results = resp(prep.out = comp_mat, model = model, weight.tgv = FALSE)
##'  }


resp = function(prep.out, model, weight.tgv = FALSE) {
  
  requireNamespace('ggplot2')
  
  prep.out <<- prep.out
  control = attr(prep.out, 'control')
  model = model
  output = list()
  fixed = model$fixed
  random = model$random
  
  ## Messages and warnings
  stopifnot("The model did not converge!" = model$converge)
  
  
  if('lrt' %in% names(model)) output$lrt = model$lrt
  
  # Variance components ---------------------
  varcomp = summary(model)$varcomp
  
  ## Dealing with the names --------------------
  if(control[,6] > 0){
    varcomp = varcomp[-grep('!R$', rownames(varcomp)),]
    rownames(varcomp)[
      rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
                                                    grepl('_1', rownames(varcomp)) &
                                                    grepl(names(control)[6],
                                                          rownames(varcomp))),])
    ] = paste("DGE", names(control)[6], sep = ':')
    rownames(varcomp)[
      rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
                                                    grepl('_2', rownames(varcomp)) &
                                                    grepl(names(control)[6],
                                                          rownames(varcomp))),])
    ] = paste("IGE", names(control)[6], sep = ':')
    
    # rownames(varcomp)[
    #   rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
    #                                                 grepl('cor', rownames(varcomp)) &
    #                                                 grepl(names(control)[6],
    #                                                       rownames(varcomp))),])
    # ] = paste("cor(IGE_DGE)", names(control)[6], sep = ':')
    
    if(control[,7] > 0){
      rownames(varcomp)[
        rownames(varcomp) %in% rownames(varcomp[which(grepl(paste(names(control)[6], 'cor$', sep='!'), 
                                                            rownames(varcomp))),])
      ] = paste(paste0("R=autocor(", 
                       paste(names(control[7]),
                             1:as.numeric(control[7]), sep = '_'), ')'),
                names(control)[6], sep = '!')
      
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl(paste(names(control)[6], 
                                             '_', sep=''),
                                       rownames(varcomp))),])
      ] = paste0("R=", rownames(varcomp[which(grepl(paste(names(control)[6], 
                                                          '_', sep=''),
                                                    rownames(varcomp))),]))
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl(paste0(names(control)[4], '!'),
                                       rownames(varcomp))),])
      ] = paste0('R=autocor(row):', 
                 sub('!.*', '', 
                     rownames(varcomp[which(grepl(names(control)[4], 
                                                  rownames(varcomp))),])))
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl(paste0(names(control)[5], "!"),
                                       rownames(varcomp))),])
      ] = paste0('R=autocor(col):', 
                 sub('!.*', '', 
                     rownames(varcomp[which(grepl(paste0(names(control)[5], "!"),
                                                  rownames(varcomp))),])))
      
    }else{
      rownames(varcomp)[
        rownames(varcomp) == rownames(varcomp[which(grepl(paste(names(control)[6], 'cor$', sep='!'), 
                                                          rownames(varcomp))),])
      ] = paste0("R=autocor(", names(control)[6],')')
      
      rownames(varcomp)[
        rownames(varcomp) %in% rownames(varcomp[which(grepl(paste(names(control)[6], '_', sep=''), 
                                                            rownames(varcomp))),])
      ] = paste0("R=", paste(names(control)[6],
                             levels(prep.out$data[,names(control)[6]]), 
                             sep = '_'))
      
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl(paste0(names(control)[4],'!'),
                                       rownames(varcomp))),])
      ] = 'R=autocor(row)'
      
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl(paste0(names(control)[5],'!'),
                                       rownames(varcomp))),])
      ] = 'R=autocor(col)'
      
    }
  }else{
    
    if(control[,7] > 0){
      rownames(varcomp)[
        rownames(varcomp) %in% rownames(varcomp[which(grepl('!R', 
                                                            rownames(varcomp))),])
      ] = paste0("R=", 
                 sub('!.*', '',
                     rownames(varcomp[which(grepl('col',rownames(varcomp))),])))
      
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('row!',
                                       rownames(varcomp))),])
      ] = paste0('R=autocor(row):', 
                 sub('!.*', '', 
                     rownames(varcomp[which(grepl('row',rownames(varcomp))),])))
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('col',
                                       rownames(varcomp))),])
      ] = paste0('R=autocor(col):', 
                 sub('!.*', '', 
                     rownames(varcomp[which(grepl('col',rownames(varcomp))),])))
      
    }else{
      rownames(varcomp)[
        rownames(varcomp) %in% rownames(varcomp[which(grepl('!R$', 
                                                            rownames(varcomp))),])
      ] = 'R'
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl(paste0(names(control)[4],'!cor'),
                                       rownames(varcomp))),])
      ] = 'R=autocor(row):'
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl(paste0(names(control)[5],'!cor'),
                                       rownames(varcomp))),])
      ] = 'R=autocor(col):'
    }
    
  }
  
  rownames(varcomp)[
    rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
                                                  grepl('_1', rownames(varcomp))),])
  ] = 'DGE'
  rownames(varcomp)[
    rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
                                                  grepl('_2', rownames(varcomp))),])
  ] = 'IGE'
  rownames(varcomp)[
    rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
                                                  grepl('cor$', rownames(varcomp))),])
  ] = "cor(IGE_DGE)"
  
  output$varcomp = varcomp
  
  
  # BLUPs --------------------
  blup = as.data.frame(summary(model, coef = TRUE)$coef.random)
  
  if(control[,6] > 0){
    ## Main effects --------------
    DGE = blup[which(grepl(names(control)[2], rownames(blup)) & 
                       !grepl(names(control)[6], rownames(blup))), -3]
    DGE[,names(control)[2]] = gsub(paste0(names(control)[2],'_'), 
                                   '', rownames(DGE), fixed = T)
    rownames(DGE) = NULL
    DGE[,names(control)[2]] = as.factor(DGE[,names(control)[2]])
    
    if(any(grepl(' ',DGE[,names(control)[2]]))){
      DGE[,names(control)[2]][grepl(' ',DGE[,names(control)[2]])] = 
        gsub(' ', '.', DGE[,names(control)[2]][grepl(' ',DGE[,names(control)[2]])])
      
    }
    
    DGE$rel.DGE = 1-(DGE$std.error^2/varcomp['DGE','component'])
    DGE = DGE[,c(3,1,2,4)]; colnames(DGE)[c(2,3)] = c('DGE', 'se.DGE')
    
    IGE = blup[which(grepl('grp', rownames(blup)) & 
                       !grepl(names(control)[6], rownames(blup))), -3]
    IGE[,names(control)[2]] = regmatches(
      do.call(rbind, strsplit(rownames(IGE), '\\(g1\\)_'))[,2], 
      m = regexpr(paste(DGE[,names(control)[2]], collapse = "|"), 
                  do.call(rbind, strsplit(rownames(IGE), '\\(g1\\)_'))[,2])
    )
    rownames(IGE) = NULL
    IGE$rel.IGE = 1-(IGE$std.error^2/varcomp['IGE','component'])
    IGE = IGE[,c(3,1,2,4)]; colnames(IGE)[c(2,3)] = c('IGE', 'se.IGE')
    
    IGE$class = ifelse(IGE$IGE > mean(IGE$IGE) + stats::sd(IGE$IGE), 'Sensitive', 
                       ifelse(IGE$IGE < mean(IGE$IGE) - stats::sd(IGE$IGE), 'Aggressive', 
                              'Homeostatic'))
    
    main = merge(DGE, IGE, by = names(control)[2])
    
    if(isFALSE(weight.tgv)){
      main$TGV = main$DGE + mean(do.call(c,lapply(Filter(function(x) is.list(x), 
                                                         prep.out), function(x){x$CIF}))) *
        main$IGE
    }else if(is.numeric(weight.tgv)){
      main$TGV = weight.tgv[1] * main$DGE + 
        mean(do.call(c,lapply(Filter(function(x) is.list(x),
                                     prep.out), function(x){x$CIF}))) *
        main$IGE * weight.tgv[2]
    }else if(isTRUE(weight.tgv)){
      main$TGV = main$rel.DGE * main$DGE + 
        mean(do.call(c,lapply(Filter(function(x) is.list(x),
                                     prep.out), function(x){x$CIF}))) *
        main$IGE * main$rel.IGE
    }
    
    ## Interaction effects ----------------
    DGE.int = blup[which(grepl(names(control)[2], rownames(blup)) & 
                           grepl(names(control)[6], rownames(blup))), -3]
    DGE.int[,names(control)[2]] = gsub(paste0(names(control)[2],'_'),'', 
                                       gsub(':.*', '', rownames(DGE.int)), fixed = T)
    DGE.int[,names(control)[6]] = gsub(paste0(names(control)[6], '_'),'', 
                                       gsub('.*:', '', rownames(DGE.int)), fixed = T)
    rownames(DGE.int) = NULL
    
    if(any(grepl(' ',DGE.int[,names(control)[2]]))){
      DGE.int[,names(control)[2]][grepl(' ',DGE.int[,names(control)[2]])] = 
        gsub(' ', '.', DGE.int[,names(control)[2]][grepl(' ',DGE.int[,names(control)[2]])])
      
    }
    
    DGE.int$rel.DGE = 1-(DGE.int$std.error^2/varcomp['DGE','component'])
    DGE.int = DGE.int[,c(3,4,1,2,5)]; colnames(DGE.int)[c(3,4)] = c('DGE', 'se.DGE')
    DGE.int = merge(DGE.int, DGE[,c(names(control)[2], 'DGE')],
                    by = names(control)[2])
    DGE.int$DGE = DGE.int$DGE.x + DGE.int$DGE.y
    DGE.int = DGE.int[,c(1, 2, 7, 4, 5)]
    
    IGE.int = blup[which(grepl('grp', rownames(blup)) & 
                           grepl(names(control)[6], rownames(blup))), -3]
    IGE.int[,names(control)[2]] = regmatches(
      do.call(rbind, strsplit(rownames(IGE.int), '\\(g1\\)_'))[,2], 
      m = regexpr(paste(DGE[,names(control)[2]], collapse = "|"), 
                  do.call(rbind, strsplit(rownames(IGE.int), '\\(g1\\)_'))[,2])
    )
    IGE.int[,names(control)[6]] = gsub(paste0(names(control)[6], '_'),'', 
                                       gsub('.*:', '', rownames(IGE.int)), fixed = T)
    rownames(IGE.int) = NULL
    IGE.int$rel.IGE = 1-(IGE.int$std.error^2/varcomp['IGE','component'])
    IGE.int = IGE.int[,c(3,4,1,2,5)]; colnames(IGE.int)[c(3,4)] = c('IGE', 'se.IGE')
    IGE.int = merge(IGE.int, IGE[,c(names(control)[2], 'IGE')],
                    by = names(control)[2])
    IGE.int$IGE = IGE.int$IGE.x + IGE.int$IGE.y
    IGE.int = IGE.int[,c(1, 2, 7, 4, 5)]
    IGE.int$class = do.call(c, lapply(split(IGE.int, IGE.int[,2]), function(x){
      ifelse(x$IGE > mean(x$IGE) + stats::sd(x$IGE), 'Sensitive', 
             ifelse(x$IGE < mean(x$IGE) - stats::sd(x$IGE), 'Aggressive', 
                    'Homeostatic'))
    }))
    
    within = merge(DGE.int, IGE.int, 
                   by = c(names(control)[2],
                          names(control)[6]))
    
    if(isFALSE(weight.tgv)){
      within =  merge(within, do.call(rbind, lapply(split(within, within[,2]), function(x){
        data.frame(
          gen = x[,1], 
          age = x[,2],
          TGV = x$DGE + prep.out[[grep(unique(x[,2]),names(prep.out))]]$CIF *
            x$IGE
        )
      })), by.x = c(names(control)[2], names(control)[6]), 
      by.y = c('gen','age'))
    }else if(is.numeric(weight.tgv)){
      within =  merge(within, do.call(rbind, lapply(split(within, within[,2]), function(x){
        data.frame(
          gen = x[,1], 
          age = x[,2],
          TGV =  weight.tgv[1] * x$DGE + 
            prep.out[[grep(unique(x[,2]),names(prep.out))]]$CIF *
            x$IGE* weight.tgv[2]
        )
      })), by.x = c(names(control)[2], names(control)[6]), 
      by.y = c('gen','age'))
    }else if(isTRUE(weight.tgv)){
      within =  merge(within, do.call(rbind, lapply(split(within, within[,2]), function(x){
        data.frame(
          gen = x[,1], 
          age = x[,2],
          TGV =  x$rel.DGE * x$DGE + 
            prep.out[[grep(unique(x[,2]),names(prep.out))]]$CIF *
            x$IGE* x$rel.IGE
        )
      })), by.x = c(names(control)[2], names(control)[6]), 
      by.y = c('gen','age'))
    }
    
    output$blups = list(main = main, within = within)
    
    ## Other random effects ----------------
    if(any(!grepl("grp\\(g1\\)", attr(model$formulae$random, 'term.labels')))) {
      coef.random = blup[which(!grepl('grp', rownames(blup)) & 
                                 !grepl(names(control)[2], rownames(blup))), -3]
      output$blups$coef.random = coef.random
    }
  }else{
    
    ## Main effects --------------
    DGE = blup[which(grepl(names(control)[2], rownames(blup))), -3]
    DGE[,names(control)[2]] = gsub(paste0(names(control)[2],'_'), 
                                   '', rownames(DGE), fixed = T)
    rownames(DGE) = NULL
    DGE[,names(control)[2]] = as.factor(DGE[,names(control)[2]])
    
    if(any(grepl(' ',DGE[,names(control)[2]]))){
      DGE[,names(control)[2]][grepl(' ',DGE[,names(control)[2]])] = 
        gsub(' ', '.', DGE[,names(control)[2]][grepl(' ',DGE[,names(control)[2]])])
      
    }
    
    DGE$rel.DGE = 1-(DGE$std.error^2/varcomp['DGE','component'])
    DGE = DGE[,c(3,1,2,4)]; colnames(DGE)[c(2,3)] = c('DGE', 'se.DGE')
    
    IGE = blup[which(grepl('grp', rownames(blup))), -3]
    IGE[,names(control)[2]] = regmatches(
      do.call(rbind, strsplit(rownames(IGE), '\\(g1\\)_'))[,2], 
      m = regexpr(paste(DGE[,names(control)[2]], collapse = "|"), 
                  do.call(rbind, strsplit(rownames(IGE), '\\(g1\\)_'))[,2])
    )
    rownames(IGE) = NULL
    IGE$rel.IGE = 1-(IGE$std.error^2/varcomp['IGE','component'])
    IGE = IGE[,c(3,1,2,4)]; colnames(IGE)[c(2,3)] = c('IGE', 'se.IGE')
    
    IGE$class = ifelse(IGE$IGE > mean(IGE$IGE) + stats::sd(IGE$IGE), 'Sensitive', 
                       ifelse(IGE$IGE < mean(IGE$IGE) - stats::sd(IGE$IGE), 'Aggressive', 
                              'Homeostatic'))
    
    main = merge(DGE, IGE, by = names(control)[2])    
    
    if(isFALSE(weight.tgv)){
      main$TGV = main$DGE + prep.out$CIF *  main$IGE
    }else if(is.numeric(weight.tgv)){
      main$TGV = weight.tgv[1] * main$DGE + prep.out$CIF *
        main$IGE * weight.tgv[2]
    }else if(isTRUE(weight.tgv)){
      main$TGV = main$rel.DGE * main$DGE + prep.out$CIF *
        main$IGE * main$rel.IGE
    }
    
    ## Other random effects ----------------
    if(any(!grepl("grp\\(g1\\)",attr(model$formulae$random, 'term.labels')))) {
      coef.random = blup[which(!grepl('grp', rownames(blup)) & 
                                 !grepl(names(control)[2], rownames(blup))), -3]
      output$blups = list(main = main, coef.random = coef.random)
    }else{
      output$blups = list(main = main)
    }
  }
  
  attr(output, 'control') = control
  attr(output, 'data') = prep.out$data
  attr(output, 'residuals') = model$residuals
  remove(prep.out, envir = .GlobalEnv)
  
  class(output) = 'comresp'
  
  return(output)
  
}




#' Plots for the `comresp` object
#'
#' Build plots using the outputs stored in the `comresp` object.
#'
#'
#' @param object An object of class `comresp`.
#' @param category A string indicating which plot to build. See options in the Details section.
#' @param level A string indicating the information level to be used for building
#' the plots. Options are `"main"` for focusing on the main effects, or `"within"` to 
#' focus on the within-age effects effects (main + interaction effects). Defaults
#' to `"main"`, which also serves when data has a single age.
#' @param age If `level = 'within'`, a string indicating if plots should be built per age. 
#' Options are `"all"` for plotting all ages, or the name of a specific age. Defaults
#' to `"all"`, which serves when data has a single age.
#' @param ... Currently not used.
#' 
#' @method plot comresp
#' 
#' @details The available options are:
##'   \itemize{
##'   \item `class`: a density plot with the IGE distribution. The area within the
##'   distribution is filled according to the competition class.
##'   \item `DGEvIGE`: a scatter plot illustrating the relation between IGE (\emph{x}-axis) 
##'   and DGE (\emph{y}-axis). The dots are coloured according to the competition class
##'   \item `DGE.IGE`: lollipop plots representing the DGE and IGE of each genotype. The plots 
##'   are in descending order according to the DGE. The dots' colour depicts the DGE's and IGE's 
##'   reliability of each genotype.
##'   \item `TGV`: lollipop plots with the TGV of each genotype, in increasing order. 
##'   \item `nneigh`: a bar plot depicting the number of different genotypes as neighbours
##'    (total and per competition class) of each selection candidate.
##'    \item `grid.res`: a heatmap representing the grid. The cells are filled according
##'    to the residual value of each plot.
##'    \item `grid.dge`: a heatmap representing the grid. The cells are filled according
##'    to the DGE value of each genotype.
##'    \item `grid.ige`: a heatmap representing the grid. The cells are filled according
##'    to the IGE value of each genotype. 
##'    \item `grid.class`: a heatmap representing the grid. The cells are filled according
##'    to the competition class of each genotype.
##'   }
#' 
#' 
#' @seealso  [ggplot2], [gencomp::resp]
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom rlang .data
#' @importFrom stats reorder sd density
#' @importFrom ggpubr ggarrange
#' 
#' @export
#' 
#' @examples
#' \donttest{
#'  library(gencomp)
#'  comp_mat = prep(data = euca, gen = 'clone', repl = 'block', area = 'area', 
#'                  ind = 'tree', age = 'age', row = 'row', col = 'col', 
#'                  dist.col = 3, dist.row = 2, trait = 'mai', method = 'SK',
#'                  n.dec = 3, verbose = TRUE)
#'  model = asr(prep.out = comp_mat, 
#'              fixed = mai~ age, 
#'              random = ~ block:age, 
#'              cor = TRUE, maxit = 20,
#'              lrtest = FALSE)
#'  results = resp(prep.out = comp_mat, model = model, weight.tgv = FALSE)
#'  
#'  plot(results, category = 'DGEvIGE', level = 'within', age = '6y')
#'  plot(results, category = 'grid.res', level = 'main', age = '3y')
#'  plot(results, category = 'nneigh', level = 'main')
#'  plot(results, category = 'DGE.IGE', level = 'within', age ='all')
#'  plot(results, category = 'grid.dge', level = 'within', age ='3y')
#'  plot(results, category = 'grid.ige', level = 'within', age ='6y')
#'  plot(results, category = 'grid.class', level = 'within', age ='all')
#'  plot(results, category = 'class', level = 'main')
#'  plot(results, category = 'TGV', level = 'within')
#'  # Note that the ages are labelled as "3y" and "6y" in the example dataset 
#'  }

plot.comresp = function(object, category = 'DGE.IGE', level = 'main', age = 'all', ...){
  
  ### Messages and warnings
  stopifnot("The object must be of class 'comresp'" = class(object) == 'comresp')
  stopifnot("'category' should be of size 1" = length(category) == 1)
  stopifnot("Please, choose between the available categories (see Details section in 'help(plot.comresp)')" = category %in% 
              c('class', 'DGEvIGE', 'DGE.IGE', 'TGV', 'nneigh', 'grid.res', 'grid.dge', 'grid.ige', 'grid.class'))
  stopifnot("'level' should be of size 1" = length(level) == 1)
  stopifnot("Please, choose between the available levels ('main' or 'within')" = level %in% c('main','within'))


  control = attr(object, 'control')
  main = object$blups$main
  dat = attr(object, 'data')
  rownames(dat) = NULL
  dat$resid = c(attr(object, 'residuals'))
  
  stopifnot("'age' should be of size 1" = length(age) == 1)
  stopifnot("'age' does not exist" = age %in% c('all', levels(dat[,colnames(control)[6]])))
  
  if(control[,6] == 0){
    age = 'all'
    level = 'main'
  }
  
  
  if(level == 'main'){
    
    if(category == 'class'){  ### Density: IGE classification 
      den.ige.main = ggplot(data = main) + 
        geom_density(aes(x = .data$IGE, after_stat(density), fill = 'Homeostatic'), 
                     alpha = .8)
      d = ggplot2::ggplot_build(den.ige.main)$data[[1]]
      den.ige.main + geom_area(data = subset(d, d$x > mean(main$IGE) + stats::sd(main$IGE)), 
                               aes(x = .data$x, y = .data$y, fill = 'Sensitive'), alpha = 0.8) + 
        geom_area(data = subset(d, subset = d$x < mean(main$IGE) - stats::sd(main$IGE)), 
                  aes(x = .data$x, y = .data$y, fill = 'Aggressive'), alpha = .8) +
        geom_density(data = main, aes(x = .data$IGE, after_stat(density)), 
                     color = 'black', linewidth = 1.2, show.legend = F) +
        scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
        labs(x = "IGE", y = 'Density', fill = 'Class') + theme_minimal() +
        theme(legend.position = 'top')
      
    } else if(category == 'DGEvIGE'){  ### DGE vs IGE
      ggplot(data = main, aes(x = .data$IGE, y = .data$DGE)) + 
        geom_point(aes(fill = .data$class), colour = 'black', pch = 21, size = 2.2) +
        scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
        labs(fill = 'Class') + theme_minimal() +
        theme(legend.position = 'top') 
      
    } else if(category == 'DGE.IGE'){   ### DGE and IGE
      temp = data.frame(
        gen = rep(main[, 1], times = 2), 
        blup = c(main$DGE, main$IGE),
        se = c(main$se.DGE, main$se.IGE),
        rel = c(main$rel.DGE, main$rel.IGE),
        comp = rep(c('DGE', 'IGE'), each = control[,2])
      )
      
      temp$gen = factor(temp$gen, levels = temp[order(temp$blup, decreasing = TRUE)[
        order(temp$blup, decreasing = TRUE) %in% which(temp$comp == 'DGE')], 1
      ])
      
      ggplot(data = temp, aes(x = .data$gen, y = .data$blup)) + 
        facet_wrap(.~.data$comp, ncol = 1, scales = 'free_y') + 
        geom_segment(aes(x = .data$gen,
                         xend = .data$gen, y = 0, yend = .data$blup),
                     linewidth = 1) + 
        geom_point(aes(fill = rel),
                   size = 2, color = 'black',
                   alpha = 0.7, shape = 21, stroke = .4) + 
        theme(axis.text.x = element_text(angle = 90, size = 8), 
              legend.position = 'top', panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey')) +
        scale_fill_viridis_c(option = 'turbo') + 
        labs(x = names(control)[2], y = 'BLUP', fill = "Reliability") + 
        guides(fill = guide_colourbar(title.position = 'top', title.hjust = .5, 
                                      barwidth = 9, frame.colour = 'black'))
      
    } else if(category == 'TGV'){   ### TGV
      temp = data.frame(gen = main[,1], tgv = main$TGV)
      ggplot(data = temp) + 
        geom_segment(aes(x = stats::reorder(.data$gen, .data$tgv),
                         xend = .data$gen, y = 0, yend = .data$tgv),
                     linewidth = 1) + 
        geom_point(aes(x = stats::reorder(.data$gen, .data$tgv), y = .data$tgv),
                   size = 2, color = 'black', alpha = 0.7, shape = 21, stroke = .4, 
                   fill = 'darkred') + 
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, size = 8)) + 
        labs(x = names(control)[2], y = 'TGV')
      
    }else if(category == 'nneigh'){  ### Number of neighbors
      
      Z = dat[, 1:control[,2]]
      Z = ifelse(crossprod(ifelse(Z == 0, 0, 1)) == 0, 0, 1)
      diag(Z) = 0
      Z = Z[order(rownames(Z)), order(colnames(Z))]
      
      neigh = do.call(rbind, lapply(
        apply(Z, 1, function(y){
          aa = data.frame(
            neigh = regmatches(
              names(y), 
              m = regexpr(paste(main[,names(control)[2]], collapse = "|"), 
                          names(y))
            )[which(y == 1)])
            
          merge(aa, main[,c(names(control)[2],'class')], by.x = 'neigh', 
                by.y = names(control)[2])
        }), function(w){
          data.frame(table(w$class))
        }
      ))
      
      neigh$gen =  regmatches(
        rownames(neigh), 
        m = regexpr(paste(main[,names(control)[2]], collapse = "|"), 
                    rownames(neigh))
      )
      
      ggplot(data = neigh, aes(x = .data$gen, y = .data$Freq, fill = .data$Var1)) +
        geom_bar(stat = "identity", color = 'black') +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, size = 8),
              axis.title = element_text(face = 'bold'),
              legend.position = 'bottom') +
        labs(x = names(control)[2], 
             y = "No. of different genotypes as neighbours", fill = 'Class') +
        scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6'))
      
    }else if(category == 'grid.res'){
      
      if(control[7] > 0){ # Multi-areas
        
        temp = dat[,c(names(control)[c(2,4,5,7,1)],'resid')]
        
        facet.label.col = paste(names(control)[7], unique(dat[,names(control)[7]]))
        names(facet.label.col) = unique(dat[,names(control)[7]])
        
        colnames(temp) = c('gen', 'row', 'col', 'area', 'y','resid')
        
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$resid)) + 
          geom_tile(color = 'black') + 
          facet_grid(cols = vars(.data$area), 
                     scales = 'free_x', 
                     labeller = labeller(.cols = facet.label.col)) +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = 'Residual') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))
        
      }else{  # Single-area
        
        temp = dat[,c(names(control)[c(2,4,5,1)],'resid')]
        
        colnames(temp) = c('gen', 'row', 'col', 'y','resid')
        
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$resid)) + 
          geom_tile(color = 'black') + 
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = 'Residual') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))
        
      }
      
    }else if(category == 'grid.dge'){
      
      if(control[7] > 0){ # Multi-areas
        
        temp = dat[,names(control)[c(2,4,5,7,1)]]
        colnames(temp) = c('gen', 'row', 'col', 'area', 'y')
        
        temp = merge(temp, main, by.x = 'gen', by.y = names(control)[2])
        
        facet.label.col = paste(names(control)[7], unique(dat[,names(control)[7]]))
        names(facet.label.col) = unique(dat[,names(control)[7]])
        
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$DGE)) + 
          geom_tile(color = 'black') + 
          facet_grid(cols = vars(.data$area), 
                     scales = 'free_x', 
                     labeller = labeller(.cols = facet.label.col)) +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = 'DGE') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))
        
      }else{  # Single-area
        
        temp = dat[,names(control)[c(2,4,5,1)]]
        colnames(temp) = c('gen', 'row', 'col', 'y')
        
        temp = merge(temp, main, by.x = 'gen', by.y = names(control)[2])
        
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$DGE)) + 
          geom_tile(color = 'black') + 
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = 'DGE') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))
        
      }
      
    }else if(category == 'grid.ige'){
      if(control[7] > 0){ # Multi-areas
        
        temp = dat[,names(control)[c(2,4,5,7,1)]]
        colnames(temp) = c('gen', 'row', 'col', 'area', 'y')
        
        temp = merge(temp, main, by.x = 'gen', by.y = names(control)[2])
        
        facet.label.col = paste(names(control)[7], unique(dat[,names(control)[7]]))
        names(facet.label.col) = unique(dat[,names(control)[7]])
        
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$IGE)) + 
          geom_tile(color = 'black') + 
          facet_grid(cols = vars(.data$area), 
                     scales = 'free_x', 
                     labeller = labeller(.cols = facet.label.col)) +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = 'IGE') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))
        
      }else{  # Single-area
        
        temp = dat[,names(control)[c(2,4,5,1)]]
        colnames(temp) = c('gen', 'row', 'col', 'y')
        
        temp = merge(temp, main, by.x = 'gen', by.y = names(control)[2])
        
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$IGE)) + 
          geom_tile(color = 'black') + 
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = 'IGE') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))
        
      }
      
    }else if(category == 'grid.class'){
      if(control[7] > 0){ # Multi-areas
        
        temp = dat[,names(control)[c(2,4,5,7,1)]]
        colnames(temp) = c('gen', 'row', 'col', 'area', 'y')
        
        temp = merge(temp, main, by.x = 'gen', by.y = names(control)[2])
        
        facet.label.col = paste(names(control)[7], unique(dat[,names(control)[7]]))
        names(facet.label.col) = unique(dat[,names(control)[7]])
        
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$class)) + 
          geom_tile(color = 'black') + 
          facet_grid(cols = vars(.data$area), 
                     scales = 'free_x', 
                     labeller = labeller(.cols = facet.label.col)) +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6'))  + 
          labs(x = "Row", y = 'Column', fill = 'Class') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))
        
      }else{  # Single-area
        
        temp = dat[,names(control)[c(2,4,5,1)]]
        colnames(temp) = c('gen', 'row', 'col', 'y')
        
        temp = merge(temp, main, by.x = 'gen', by.y = names(control)[2])
        
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$class)) + 
          geom_tile(color = 'black') +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6'))  + 
          labs(x = "Row", y = 'Column', fill = 'Class') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))
        
      }
      
    }
    
  } else if(level == 'within'){
    within = object$blups$within
    
    if(age == 'all'){
      temp = split(within, within[,2])
      if(category == 'class'){
        
        do.call(
          ggpubr::ggarrange, 
          c(lapply(temp, function(df){
            aa = ggplot(data = df) + 
              geom_density(aes(x = .data$IGE, after_stat(density), fill = 'Homeostatic'), 
                           alpha = .8)
            d = ggplot2::ggplot_build(aa)$data[[1]]
            aa = aa + geom_area(data = subset(d, d$x > mean(df$IGE) + stats::sd(df$IGE)), 
                                aes(x = .data$x, y = .data$y, fill = 'Sensitive'), alpha = 0.8) + 
              geom_area(data = subset(d, d$x < mean(df$IGE) - stats::sd(df$IGE)), 
                        aes(x = .data$x, y = .data$y, fill = 'Aggressive'), alpha = .8) +
              geom_density(data = df, aes(x = .data$IGE, after_stat(density)), 
                           color = 'black', linewidth = 1.2, show.legend = F) +
              scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
              labs(x = "IGE", y = 'Density', fill = 'Class', 
                   subtitle = paste(names(control)[6], unique(df[,2]))) + 
              theme_minimal() +
              theme(legend.position = 'top')
          }), common.legend = TRUE)
        )
        
      }else if(category == 'DGEvIGE'){
        
        facet.label =  paste(names(control)[6], unique(within[,names(control)[6]]))
        names(facet.label) = unique(within[,names(control)[6]])
        
        ggplot(data = within, aes(x = .data$IGE, y = .data$DGE)) + 
          geom_point(aes(fill = .data$class), colour = 'black', pch = 21, size = 2.2) +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
          labs(fill = 'Class') +
          theme_bw() +
          theme(legend.position = 'top') +
          facet_wrap(.~.data$age, scales = 'free', labeller = labeller(.cols = facet.label))
        
        
      }else if(category == 'DGE.IGE'){
        
        do.call(
          ggpubr::ggarrange, 
          c(lapply(temp, function(df){
            
            temp2 = data.frame(
              gen = rep(df[, 1], times = 2), 
              blup = c(df$DGE, df$IGE),
              se = c(df$se.DGE, df$se.IGE),
              rel = c(df$rel.DGE, df$rel.IGE),
              comp = rep(c('DGE', 'IGE'), each = control[,2])
            )
            
            temp2$gen = factor(temp2$gen, levels = temp2[order(temp2$blup, decreasing = TRUE)[
              order(temp2$blup, decreasing = TRUE) %in% which(temp2$comp == 'DGE')], 1
            ])
            
            aa = ggplot(data = temp2, aes(x = .data$gen, y = .data$blup)) + 
              facet_wrap(.~.data$comp, ncol = 1, scales = 'free_y') + 
              geom_segment(aes(x = .data$gen,
                               xend = .data$gen, y = 0, yend = .data$blup),
                           linewidth = 1) + 
              geom_point(aes(fill = rel),
                         size = 2, color = 'black',
                         alpha = 0.7, shape = 21, stroke = .4) + 
              theme_bw()+
              theme(axis.text.x = element_text(angle = 90, size = 8), 
                    legend.position = 'top') +
              scale_fill_viridis_c(option = 'turbo') + 
              labs(x = names(control)[2], y = 'BLUP', fill = "Reliability") + 
              guides(fill = guide_colourbar(title.position = 'top', title.hjust = .5, 
                                            barwidth = 9, frame.colour = 'black'))
            
          }), common.legend = TRUE)
        )
        
      }else if(category == 'TGV'){
        
        temp3 = data.frame(gen = within[,1], tgv = within$TGV, 
                           age = within[,2])
        
        retrieve = function(x) do.call(rbind, strsplit(x, '@#_'))[,1]
        
        facet.label =  paste(names(control)[6], unique(within[,names(control)[6]]))
        names(facet.label) = unique(within[,names(control)[6]])
        
        
        ggplot(data = cbind(temp3, V4 = paste(temp3$gen, temp3$age, sep = '@#_')),
               aes(x = stats::reorder(.data$V4, -.data$tgv), y = .data$tgv)) +
          facet_wrap(.~.data$age, scales = 'free', labeller = labeller(.cols = facet.label)) + 
          theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 8))+
          scale_x_discrete(labels = retrieve) + 
          geom_segment(aes(x = stats::reorder(.data$V4, .data$tgv),
                           xend = .data$V4, y = 0, yend = .data$tgv),
                       linewidth = 1) + 
          geom_point(aes(x = stats::reorder(.data$V4, .data$tgv), y = .data$tgv),
                     size = 2, color = 'black', alpha = 0.7, shape = 21, stroke = .4, 
                     fill = 'darkred') + 
          labs(x = names(control)[2], y = 'TGV') 
        
        
      }else if(category == 'nneigh'){
        
        temp = split(dat, dat[,names(control)[6]])
        
        
        neigh = do.call(rbind,lapply(temp, function(x){
          Z = x[, 1:control[,2]]
          Z = ifelse(crossprod(ifelse(Z == 0, 0, 1)) == 0, 0, 1)
          diag(Z) = 0
          Z = Z[order(rownames(Z)), order(colnames(Z))]
          
          neigh = do.call(rbind, lapply(
            apply(Z, 1, function(y){
              data.frame(
                neigh = regmatches(
                  names(y), 
                  m = regexpr(paste(within[,names(control)[2]], collapse = "|"), 
                              names(y))
                )[which(y == 1)],
                class = within[
                  which(within[,2] == unique(x[, names(control)[6]])),
                  'class'
                ][
                  match(regmatches(
                    names(y), 
                    m = regexpr(paste(within[,names(control)[2]], collapse = "|"), 
                                names(y))
                  )[which(y == 1)],
                  within[which(within[,2] == unique(x[, names(control)[6]])),
                         c(names(control)[2])])
                ]
              )
            }), function(w){
              data.frame(table(w$class))
            }
          ))
          
          neigh$gen =  regmatches(
            rownames(neigh), 
            m = regexpr(paste(within[,names(control)[2]], collapse = "|"), 
                        rownames(neigh))
          )
          neigh$age = unique(x[,names(control)[6]])
          neigh
        }))
        
        facet.label =  paste(names(control)[6], unique(within[,names(control)[6]]))
        names(facet.label) = unique(within[,names(control)[6]])
        
        ggplot(data = neigh, aes(x = .data$gen, y = .data$Freq, fill = .data$Var1)) +
          geom_bar(stat = "identity", color = 'black') +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, size = 8),
                axis.title = element_text(face = 'bold'),
                legend.position = 'bottom') +
          labs(x = names(control)[2], 
               y = "No. of different genotypes as neighbours", fill = 'Class') +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
          facet_wrap(.~.data$age, labeller = labeller(.cols = facet.label))
      }else if(category == 'grid.res'){
        
        if(control[7] > 0){ # Multi-areas
          
          temp2 = dat[,c(names(control)[c(2,4,5,6,7,1)],'resid')]
          
          facet.label.col = paste(names(control)[7], unique(dat[,names(control)[7]]))
          names(facet.label.col) = unique(dat[,names(control)[7]])
          facet.label.row =  paste(names(control)[6], unique(dat[,names(control)[6]]))
          names(facet.label.row) = unique(dat[,names(control)[6]])
          
          colnames(temp2) = c('gen', 'row', 'col', 'age', 'area', 'y','resid')
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$resid)) + 
            geom_tile(color = 'black') + 
            facet_grid(cols = vars(.data$area), rows = vars(.data$age),
                       scales = 'free_x', 
                       labeller = labeller(.cols = facet.label.col,
                                           .rows = facet.label.row)) +
            scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
            labs(x = "Row", y = 'Column', fill = 'Residual') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }else{  # Single-area
          
          temp2 = dat[,c(names(control)[c(2,4,5,6,1)],'resid')]
          
          facet.label.row =  paste(names(control)[6], unique(dat[,names(control)[6]]))
          names(facet.label.row) = unique(dat[,names(control)[6]])
          
          colnames(temp2) = c('gen', 'row', 'col', 'age', 'y','resid')
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$resid)) + 
            geom_tile(color = 'black') + 
            facet_grid(rows = vars(.data$age),
                       scales = 'free_x', 
                       labeller = labeller(.rows = facet.label.row)) +
            scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
            labs(x = "Row", y = 'Column', fill = 'Residual') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }
        
      }else if(category == 'grid.dge'){
        
        if(control[7] > 0){ # Multi-area
          
          temp2 = dat[,names(control)[c(2,4,5,6,7,1)]]
          colnames(temp2) = c('gen', 'row', 'col', 'age', 'area', 'y')
          
          temp2 = merge(temp2, within, by.x = c('gen','age'), 
                        by.y = c(names(control)[2], names(control)[6]))
          
          facet.label.col = paste(names(control)[7], unique(dat[,names(control)[7]]))
          names(facet.label.col) = unique(dat[,names(control)[7]])
          facet.label.row =  paste(names(control)[6], unique(dat[,names(control)[6]]))
          names(facet.label.row) = unique(dat[,names(control)[6]])
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$DGE)) + 
            geom_tile(color = 'black') + 
            facet_grid(cols = vars(.data$area), rows = vars(.data$age),
                       scales = 'free_x', 
                       labeller = labeller(.cols = facet.label.col,
                                           .rows = facet.label.row)) +
            scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
            labs(x = "Row", y = 'Column', fill = 'DGE') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }else{  # Single-area
          
          temp2 = dat[,names(control)[c(2,4,5,6,1)]]
          colnames(temp2) = c('gen', 'row', 'col', 'age', 'y')
          
          temp2 = merge(temp2, within, by.x = c('gen','age'), 
                        by.y = c(names(control)[2], names(control)[6]))
          
          facet.label.row =  paste(names(control)[6], unique(dat[,names(control)[6]]))
          names(facet.label.row) = unique(dat[,names(control)[6]])
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$DGE)) + 
            geom_tile(color = 'black') + 
            facet_grid(rows = vars(.data$age),
                       scales = 'free_x', 
                       labeller = labeller(.rows = facet.label.row)) +
            scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
            labs(x = "Row", y = 'Column', fill = 'DGE') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }
        
      }else if(category == 'grid.ige'){
        if(control[7] > 0){ # Multi-area
          
          temp2 = dat[,names(control)[c(2,4,5,6,7,1)]]
          colnames(temp2) = c('gen', 'row', 'col', 'age', 'area', 'y')
          
          temp2 = merge(temp2, within, by.x = c('gen','age'), 
                        by.y = c(names(control)[2], names(control)[6]))
          
          facet.label.col = paste(names(control)[7], unique(dat[,names(control)[7]]))
          names(facet.label.col) = unique(dat[,names(control)[7]])
          facet.label.row =  paste(names(control)[6], unique(dat[,names(control)[6]]))
          names(facet.label.row) = unique(dat[,names(control)[6]])
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$IGE)) + 
            geom_tile(color = 'black') + 
            facet_grid(cols = vars(.data$area), rows = vars(.data$age),
                       scales = 'free_x', 
                       labeller = labeller(.cols = facet.label.col,
                                           .rows = facet.label.row)) +
            scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
            labs(x = "Row", y = 'Column', fill = 'IGE') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }else{  # Single-area
          
          temp2 = dat[,names(control)[c(2,4,5,6,1)]]
          colnames(temp2) = c('gen', 'row', 'col', 'age', 'y')
          
          temp2 = merge(temp2, within, by.x = c('gen','age'), 
                        by.y = c(names(control)[2], names(control)[6]))
          
          facet.label.row =  paste(names(control)[6], unique(dat[,names(control)[6]]))
          names(facet.label.row) = unique(dat[,names(control)[6]])
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$IGE)) + 
            geom_tile(color = 'black') + 
            facet_grid(rows = vars(.data$age),
                       scales = 'free_x', 
                       labeller = labeller(.rows = facet.label.row)) +
            scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
            labs(x = "Row", y = 'Column', fill = 'IGE') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }
        
      }else if(category == 'grid.class'){
        if(control[7] > 0){ # Multi-area
          
          temp2 = dat[,names(control)[c(2,4,5,6,7,1)]]
          colnames(temp2) = c('gen', 'row', 'col', 'age', 'area', 'y')
          
          temp2 = merge(temp2, within, by.x = c('gen','age'), 
                        by.y = c(names(control)[2], names(control)[6]))
          
          facet.label.col = paste(names(control)[7], unique(dat[,names(control)[7]]))
          names(facet.label.col) = unique(dat[,names(control)[7]])
          facet.label.row =  paste(names(control)[6], unique(dat[,names(control)[6]]))
          names(facet.label.row) = unique(dat[,names(control)[6]])
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$class)) + 
            geom_tile(color = 'black') + 
            facet_grid(cols = vars(.data$area), rows = vars(.data$age),
                       scales = 'free_x', 
                       labeller = labeller(.cols = facet.label.col,
                                           .rows = facet.label.row)) +
            scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6'))  + 
            labs(x = "Row", y = 'Column', fill = 'Class') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }else{  # Single-area
          
          temp2 = dat[,names(control)[c(2,4,5,6,1)]]
          colnames(temp2) = c('gen', 'row', 'col', 'age', 'y')
          
          temp2 = merge(temp2, within, by.x = c('gen','age'), 
                        by.y = c(names(control)[2], names(control)[6]))
          
          facet.label.row =  paste(names(control)[6], unique(dat[,names(control)[6]]))
          names(facet.label.row) = unique(dat[,names(control)[6]])
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$class)) + 
            geom_tile(color = 'black') + 
            facet_grid(rows = vars(.data$age),
                       scales = 'free_x', 
                       labeller = labeller(.rows = facet.label.row)) +
            scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6'))  + 
            labs(x = "Row", y = 'Column', fill = 'Class') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }
        
      }
      
    } else {  ### Each age particularly
      
      temp = within[which(within[,names(control)[6]] == age),]
      
      
      if(category == 'class'){
        
        aa = ggplot(data = temp) + 
          geom_density(aes(x = .data$IGE, after_stat(density), fill = 'Homeostatic'), 
                       alpha = .8)
        d = ggplot2::ggplot_build(aa)$data[[1]]
        aa + geom_area(data = subset(d, d$x > mean(temp$IGE) + stats::sd(temp$IGE)), 
                       aes(x = .data$x, y = .data$y, fill = 'Sensitive'), alpha = 0.8) + 
          geom_area(data = subset(d, d$x < mean(temp$IGE) - stats::sd(temp$IGE)), 
                    aes(x = .data$x, y = .data$y, fill = 'Aggressive'), alpha = .8) +
          geom_density(data = temp, aes(x = .data$IGE, after_stat(density)), 
                       color = 'black', linewidth = 1.2, show.legend = F) +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
          labs(x = "IGE", y = 'Density', fill = 'Class') + 
          theme_minimal() +
          theme(legend.position = 'top')
        
      }else if(category == 'DGEvIGE'){
        
        ggplot(data = temp, aes(x = .data$IGE, y = .data$DGE)) + 
          geom_point(aes(fill = .data$class), colour = 'black', pch = 21, size = 2.2) +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
          labs(fill = 'Class') +
          theme_bw() +
          theme(legend.position = 'top') 
        
        
      }else if(category == 'DGE.IGE'){
        
        temp2 = data.frame(
          gen = rep(temp[, 1], times = 2), 
          blup = c(temp$DGE, temp$IGE),
          se = c(temp$se.DGE, temp$se.IGE),
          rel = c(temp$rel.DGE, temp$rel.IGE),
          comp = rep(c('DGE', 'IGE'), each = control[,2])
        )
        
        temp2$gen = factor(temp2$gen, levels = temp2[order(temp2$blup, decreasing = TRUE)[
          order(temp2$blup, decreasing = TRUE) %in% which(temp2$comp == 'DGE')], 1
        ])
        
        ggplot(data = temp2, aes(x = .data$gen, y = .data$blup)) + 
          facet_wrap(.~.data$comp, ncol = 1, scales = 'free_y') + 
          geom_segment(aes(x = .data$gen,
                           xend = .data$gen, y = 0, yend = .data$blup),
                       linewidth = 1) + 
          geom_point(aes(fill = rel),
                     size = 2, color = 'black',
                     alpha = 0.7, shape = 21, stroke = .4) + 
          theme_bw()+
          theme(axis.text.x = element_text(angle = 90, size = 8), 
                legend.position = 'top') +
          scale_fill_viridis_c(option = 'turbo') + 
          labs(x = names(control)[2], y = 'BLUP', fill = "Reliability") + 
          guides(fill = guide_colourbar(title.position = 'top', title.hjust = .5, 
                                        barwidth = 9, frame.colour = 'black'))
        
        
      }else if(category == 'TGV'){
        
        temp3 = data.frame(gen = temp[,1], tgv = temp$TGV, 
                           age = temp[,2])
        
        ggplot(data = temp3) + 
          geom_segment(aes(x = stats::reorder(.data$gen, .data$tgv),
                           xend = .data$gen, y = 0, yend = .data$tgv),
                       linewidth = 1) + 
          geom_point(aes(x = stats::reorder(.data$gen, .data$tgv), y = .data$tgv),
                     size = 2, color = 'black', alpha = 0.7, shape = 21, stroke = .4, 
                     fill = 'darkred') + 
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, size = 8)) + 
          labs(x = names(control)[2], y = 'TGV')
        
      }else if(category == 'nneigh'){  ### Number of neighbors
        
        temp = dat[which(dat[,names(control)[6]] == age),]
        
        Z = temp[, 1:control[,2]]
        Z = ifelse(crossprod(ifelse(Z == 0, 0, 1)) == 0, 0, 1)
        diag(Z) = 0
        Z = Z[order(rownames(Z)), order(colnames(Z))]
        
        neigh = do.call(rbind, lapply(
          apply(Z, 1, function(y){
            data.frame(
              neigh = regmatches(
                names(y), 
                m = regexpr(paste(within[,names(control)[2]], collapse = "|"), 
                            names(y))
              )[which(y == 1)],
              class = within[
                which(within[,2] == unique(temp[, names(control)[6]])),
                'class'
              ][
                match(regmatches(
                  names(y), 
                  m = regexpr(paste(within[,names(control)[2]], collapse = "|"), 
                              names(y))
                )[which(y == 1)],
                within[which(within[,2] == unique(temp[, names(control)[6]])),
                       c(names(control)[2])])
              ]
            )
          }), function(w){
            data.frame(table(w$class))
          }
        ))
        
        neigh$gen =  regmatches(
          rownames(neigh), 
          m = regexpr(paste(within[,names(control)[2]], collapse = "|"), 
                      rownames(neigh))
        )
        
        ggplot(data = neigh, aes(x = .data$gen, y = .data$Freq, fill = .data$Var1)) +
          geom_bar(stat = "identity", color = 'black') +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, size = 8),
                axis.title = element_text(face = 'bold'),
                legend.position = 'bottom') +
          labs(x = names(control)[2], 
               y = "No. of different genotypes as neighbours", fill = 'Class') +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6'))
        
      }else if(category == 'grid.res'){
        
        if(control[7] > 0){ # Multi-areas
          
          temp = dat[which(dat[,names(control)[6]] == age),]
          temp2 = temp[,c(names(control)[c(2,4,5,7,1)],'resid')]
          
          facet.label.col = paste(names(control)[7], unique(temp2[,names(control)[7]]))
          names(facet.label.col) = unique(temp2[,names(control)[7]])
          
          colnames(temp2) = c('gen', 'row', 'col', 'area', 'y','resid')
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$resid)) + 
            geom_tile(color = 'black') + 
            facet_grid(cols = vars(.data$area), 
                       scales = 'free_x', 
                       labeller = labeller(.cols = facet.label.col)) +
            scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
            labs(x = "Row", y = 'Column', fill = 'Residual') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }else{  # Single-area
          
          temp = dat[which(dat[,names(control)[6]] == age),]
          temp2 = dat[,c(names(control)[c(2,4,5,1)],'resid')]
          
          colnames(temp2) = c('gen', 'row', 'col', 'y','resid')
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$resid)) + 
            geom_tile(color = 'black') + 
            scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
            labs(x = "Row", y = 'Column', fill = 'Residual') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }
        
      }else if(category == 'grid.dge'){
        
        if(control[7] > 0){ # Multi-area
          
          temp3 = dat[which(dat[,names(control)[6]] == age),]
          temp2 = temp3[,names(control)[c(2,4,5,7,1)]]
          colnames(temp2) = c('gen', 'row', 'col', 'area', 'y')
          
          temp2 = merge(temp2, temp, by.x = 'gen', by.y = names(control)[2])
          
          facet.label.col = paste(names(control)[7], unique(dat[,names(control)[7]]))
          names(facet.label.col) = unique(dat[,names(control)[7]])
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$DGE)) + 
            geom_tile(color = 'black') + 
            facet_grid(cols = vars(.data$area), 
                       scales = 'free_x', 
                       labeller = labeller(.cols = facet.label.col)) +
            scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
            labs(x = "Row", y = 'Column', fill = 'DGE') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }else{  # Single-area
          
          temp3 = dat[which(dat[,names(control)[6]] == age),]
          temp2 = temp3[,names(control)[c(2,4,5,1)]]
          colnames(temp2) = c('gen', 'row', 'col', 'y')
          
          temp2 = merge(temp2, temp, by.x = 'gen', by.y = names(control)[2])
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$DGE)) + 
            geom_tile(color = 'black') + 
            scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
            labs(x = "Row", y = 'Column', fill = 'DGE') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }
        
      }else if(category == 'grid.ige'){
        if(control[7] > 0){ # Multi-areas
          
          temp3 = dat[which(dat[,names(control)[6]] == age),]
          temp2 = temp3[,names(control)[c(2,4,5,7,1)]]
          colnames(temp2) = c('gen', 'row', 'col', 'area', 'y')
          
          temp2 = merge(temp2, temp, by.x = 'gen', by.y = names(control)[2])
          
          facet.label.col = paste(names(control)[7], unique(dat[,names(control)[7]]))
          names(facet.label.col) = unique(dat[,names(control)[7]])
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$IGE)) + 
            geom_tile(color = 'black') + 
            facet_grid(cols = vars(.data$area), 
                       scales = 'free_x', 
                       labeller = labeller(.cols = facet.label.col)) +
            scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
            labs(x = "Row", y = 'Column', fill = 'IGE') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }else{  # Single-area
          
          temp3 = dat[which(dat[,names(control)[6]] == age),]
          temp2 = temp3[,names(control)[c(2,4,5,1)]]
          colnames(temp2) = c('gen', 'row', 'col', 'y')
          
          temp2 = merge(temp2, temp, by.x = 'gen', by.y = names(control)[2])
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$IGE)) + 
            geom_tile(color = 'black') + 
            scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
            labs(x = "Row", y = 'Column', fill = 'IGE') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }
        
      }else if(category == 'grid.class'){
        if(control[7] > 0){ # Multi-areas
          
          temp3 = dat[which(dat[,names(control)[6]] == age),]
          temp2 = temp3[,names(control)[c(2,4,5,7,1)]]
          colnames(temp2) = c('gen', 'row', 'col', 'area', 'y')
          
          temp2 = merge(temp2, temp, by.x = 'gen', by.y = names(control)[2])
          
          facet.label.col = paste(names(control)[7], unique(dat[,names(control)[7]]))
          names(facet.label.col) = unique(dat[,names(control)[7]])
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$class)) + 
            geom_tile(color = 'black') + 
            facet_grid(cols = vars(.data$area), 
                       scales = 'free_x', 
                       labeller = labeller(.cols = facet.label.col)) +
            scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6'))  + 
            labs(x = "Row", y = 'Column', fill = 'Class') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }else{  # Single-area
          
          temp3 = dat[which(dat[,names(control)[6]] == age),]
          temp2 = temp3[,names(control)[c(2,4,5,1)]]
          colnames(temp2) = c('gen', 'row', 'col', 'y')
          
          temp2 = merge(temp2, temp, by.x = 'gen', by.y = names(control)[2])
          
          ggplot(data = temp2, 
                 aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                     fill = .data$class)) + 
            geom_tile(color = 'black') +
            scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6'))  + 
            labs(x = "Row", y = 'Column', fill = 'Class') + 
            theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                  panel.background = element_blank(), 
                  panel.grid = element_line(colour = 'lightgrey'))
          
        }
        
      }
      
    }
    
  }
  
}


#' Print an object of class `comresp`
#'
#' Print a `comresp` object in the R console.
#'
#' @param object An object of class `comresp`.
#' @param category A string indicating which object to print. Options are `"all"`
#' for printing all objects, `"summar"` for printing the variance components and the 
#' likelihood ratio test (if `lrt = TRUE` in [gencomp::asr]), `"blup.main"` (Default) for printing
#' the DGE, IGE and TGV (the main effects if a multi-age model was fitted), and 
#' `"blup.within"` for printing the DGE, IGE and TGV within ages (if a multi-age model 
#' was fitted).
#' @param age If `category = 'blup.within'`, a string indicating if this object should be printed per age. 
#' Options are `"all"` for printing all ages, or the name of a specific age. Defaults
#' to `"all"`.
#' @param ... Currently not used.
#' 
#' @method print comresp
#' 
#' @seealso  [gencomp::resp]
#' 
#' @importFrom data.table data.table 
#' 
#' @export
#' 
#' @examples
#' \donttest{
#'  library(gencomp)
#'  comp_mat = prep(data = euca, gen = 'clone', repl = 'block', area = 'area', 
#'                  ind = 'tree', age = 'age', row = 'row', col = 'col', 
#'                  dist.col = 3, dist.row = 2, trait = 'mai', method = 'SK',
#'                  n.dec = 3, verbose = TRUE)
#'  
#'  model = asr(prep.out = comp_mat, 
#'              fixed = mai~ age, 
#'              random = ~ block:age, 
#'              cor = TRUE, maxit = 20,
#'              lrtest = FALSE)
#'  results = resp(prep.out = comp_mat, model = model, weight.tgv = FALSE)
#'  
#'  print(results)
#'  print(results, category = 'summar')
#'  print(results, category = 'blup.main')
#'  print(results, category = 'blup.within')
#'  print(results, category = 'blup.within', age = '6y')
#'  # Note that the ages are labelled as "3y" and "6y" in the example dataset 
#'  }

print.comresp = function(object, category = 'blup.main', age = 'all', ...){
  
  stopifnot("The object must be of class 'comresp'" = class(object) == 'comresp')
  stopifnot("'category' should be of size 1" = length(category) == 1)
  stopifnot("Please, choose between the available categories ('all', 'summar', 'blup.main' or 'blup.within')" = category %in% 
              c('all', 'summar', 'blup.main', 'blup.within'))

  control = attr(object, 'control')
  main = object$blups$main
  dat = attr(object, 'data')
  rownames(dat) = NULL
  dat$resid = c(attr(object, 'residuals'))
  
  stopifnot("'age' should be of size 1" = length(age) == 1)
  stopifnot("'age' does not exist" = age %in% c('all', levels(dat[,colnames(control)[6]])))
  
  if(control[,6] == 0) age = 'all'
  
  
  if(category == 'all' | category == 'summar'){

    cat("\n","===> Likelihood ratio tests", "\n")
    if('lrt' %in% names(object)) print(object$lrt)
    
    cat("\n","===> Variance components", "\n")
    print(object$varcomp)
  }
  
  
  if(category == 'all' | category == 'blup.main'){
    
    cat("\n","===> DGE, IGE and TGV", "\n")
    
    print(data.table(object$blups$main))
  }
  
  
  if(category == 'all' | category == 'blup.within'){
    
    if(age == 'all'){
      
      cat("\n","===> DGE, IGE and TGV (all ages)", "\n")
      
      print(data.table(object$blups$within))
      
      
    }else{
      
      cat("\n","===> DGE, IGE and TGV", paste0('(',age,')'), "\n")
      
      print(data.table(object$blups$within[which(object$blups$within[,2] == age),]))
      
    }
    
  }
  
}


#' Summary of the `comresp` object
#'
#' A brief summary of the `comresp` object.
#'
#' @param object An object of class `comresp`
#' @param ... Currently not used.
#' 
#' @method summary comresp
#' 
#' @return The function returns a list with the variance components (`varcomp`)
#' and a dataframe informing which genotypes had the highest DGE, the highest and 
#' the lowest IGE, and the highest TGV. This dataframe also informs the number of 
#' aggressive, homeostatic and sensitive. If a multi-age model was fitted, the output
#' list will have another dataframe with the same information for each age.
#' 
#' @seealso [gencomp::resp]
#' 
#' @export
#' 
#' @examples
#'\donttest{
#' library(gencomp)
#' comp_mat = prep(data = euca, gen = 'clone', repl = 'block', area = 'area', 
#'                 ind = 'tree', age = 'age', row = 'row', col = 'col', 
#'                 dist.col = 3, dist.row = 2, trait = 'mai', method = 'SK',
#'                 n.dec = 3, verbose = TRUE)
#'  model = asr(prep.out = comp_mat, 
#'              fixed = mai~ age, 
#'              random = ~ block:age, 
#'              cor = TRUE, maxit = 20,
#'              lrtest = FALSE)
#'  results = resp(prep.out = comp_mat, model = model, weight.tgv = FALSE)
#'  
#'  summary(results)
#' }
#'

summary.comresp = function(object, ...){
  
  output = list()
  main = object$blups$main
  
  output$varcomp = object$varcomp
  
  output$blup = data.frame(
    topDGE = as.matrix(main[order(main$DGE, decreasing = TRUE),1][1]),
    topIGE = as.matrix(main[order(main$IGE, decreasing = TRUE),1][1]),
    bottomIGE = as.matrix(main[order(main$DGE, decreasing = FALSE),1][1]),
    topTGV = main[order(main$TGV, decreasing = TRUE),'clone'][1],
    no.Aggressive = table(main$class)['Aggressive'],
    no.Homeostatic = table(main$class)['Homeostatic'],
    no.Sensitive = table(main$class)['Sensitive'],
    row.names = NULL
  )
  
  if('within' %in% names(object$blups)){
    
    within = object$blups$within
    
    temp = split(within, within[,2])
    
    temp = do.call(rbind, lapply(temp, function(df){
      
      data.frame(
        age = unique(df[,2]),
        topDGE = as.matrix(df[order(df$DGE, decreasing = TRUE),1][1]),
        topIGE = as.matrix(df[order(df$IGE, decreasing = TRUE),1][1]),
        bottomIGE = as.matrix(df[order(df$DGE, decreasing = FALSE),1][1]),
        topTGV = df[order(df$TGV, decreasing = TRUE),'clone'][1],
        no.Aggressive = table(df$class)['Aggressive'],
        no.Homeostatic = table(df$class)['Homeostatic'],
        no.Sensitive = table(df$class)['Sensitive'],
        row.names = NULL
      )
      
    }))
    rownames(temp) = NULL
    
    output$blup.within = temp
    
  }
  
return(output)
}
