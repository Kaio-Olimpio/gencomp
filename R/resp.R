##' @title Extract outputs from the genetic-spatial competition model
##' 
##' @description
##' This function provides responses extracted from the genetic-spatial 
##' competition model fitted using [competition::comp.asr()]. 
##' 
##' @param prep.out A `comp.prep` object.
##' @param model An `asreml` object, preferably obtained using the [competition::comp.asr()] function.
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
##' \item \code{blups} : The direct (DGE) and indirect genetic effects (IGE), their standard errors, 
##' the competition class of each genotype and the total genotypic value (TGV). If `age = TRUE`, 
##' `blup` is a list containing the main effects and the within-ages DGE and IGE. If 
##' other random effects were declared in the model, `blup` will contain another data frame 
##' with their BLUPs.
##' \item \code{plots} : A list of plots:
##'   \itemize{
##'   \item `IGE.density`: a density plot with the IGE distribution. The area within the
##'   distribution is filled according to the competition class (see Details).
##'   \item `DGEvsIGE`: a scatter plot illustrating the relation between IGE (\emph{x}-axis) 
##'   and DGE (\emph{y}-axis). The dots are coloured according to the competition class
##'   \item `DGE.IGE`: a lollipop plots representing the DGE and IGE of each genotype. The plots 
##'   are in descending order according to the DGE. The dots' colour depicts the DGE's and IGE's 
##'   reliability of each genotype.
##'   \item `TGV`: lollipop plots with the TGV of each genotype, in increasing order. 
##'   \item `n.neigh`: a bar plot depicting the number of different genotypes as neighbours
##'    (total and per competition class) of each selection candidate. If `age = TRUE`, 
##'    there will be a `n.neigh` plot for each age. 
##'    \item `grid.pheno`: a heatmap representing the grid. The cells are filled according
##'    to the phenotype value of each plot. If `age = TRUE`, there will be a `grid.pheno` plot for each age.
##'    \item `grid.DGE`: a heatmap representing the grid. The cells are filled according
##'    to the DGE value of each genotype. If `age = TRUE`, there will be a `grid.DGE` plot for each age.
##'    \item `grid.IGE`: a heatmap representing the grid. The cells are filled according
##'    to the IGE value of each genotype. If `age = TRUE`, there will be a `grid.DGE` plot for each age.
##'    \item `grid.IGE.class`: a heatmap representing the grid. The cells are filled according
##'    to the competition class of each genotype. If `age = TRUE`, there will be a `grid.DGE` plot for each age. 
##'   }
##'   All plots are customizable using resources of the `ggplot2` library.
##' }
##'
##' @details
##' The genetic-spatial competition model provides the direct (DGE) and indirect 
##' genetic effects (IGE) of each genotype. The DGE represents the "pure" performance
##' of the genotype, while the IGE is the related to the average effect of the genotype on
##' the genotypic value of its neighbours. The higher the IGE, the more aggressive is the genotype. 
##' Here, we use the classification proposed by \insertCite{ferreira_novel_2023;textual}{competition} 
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
##' [competition::comp.prep()]. If `weight.tgv = TRUE`, the DGE and IGE will be 
##' multiplied by their respective reliabilities (\eqn{r_{d_i}^2} and \eqn{r_{c_i}^2}):
##' 
##' \deqn{wTGV_i = r_{d_i}^2 \times d_i + r_{c_i}^2 \times {CIF \times c_i}}
##'
##' @seealso  [competition::comp.prep], [competition::comp.asr], [ggplot2]
##' 
##' @references 
##' \insertAllCited{}
##' 
##' @import ggplot2
##' @importFrom rlang .data
##' @importFrom stats reorder sd density
##' @importFrom ggpubr ggarrange
##' @importFrom Rdpack reprompt
##' @export
##' 
##' @examples
##' \donttest{
##'  comp_mat = comp.prep(data = euca, gen = 'clone', repl = 'block', area = 'area', 
##'                       ind = 'tree', age = 'age', row = 'row', col = 'col', 
##'                       dist.col = 3, dist.row = 2, trait = 'mai', method = 'SK',
##'                       n.dec = 3, verbose = TRUE)
##'  
##'  model = comp.asr(prep.out = comp_mat, 
##'                   fixed = mai~ age, 
##'                   random = ~ block:age, 
##'                   cor = TRUE, maxit = 50)
##'                   
##'   results = comp.resp(prep.out = comp_mat, model = model, weight.tgv = FALSE)
##'  }


comp.resp = function(prep.out, model, weight.tgv = FALSE) {
  
  requireNamespace('ggplot2')
  
  prep.out <<- prep.out
  model = model
  output = list()
  fixed = model$fixed
  random = model$random

  if(!model$converge) stop('The model did not converge!')
  
  output$lrt = model$lrt
  
  # Variance components ---------------------
  varcomp = summary(model)$varcomp
  
  ## Dealing with the names --------------------
  if(prep.out$control[,6] > 0){
    varcomp = varcomp[-grep('!R', rownames(varcomp)),]
    rownames(varcomp)[
      rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
                                                    grepl('_1', rownames(varcomp)) &
                                                    grepl(names(prep.out$control)[6],
                                                          rownames(varcomp))),])
    ] = paste("DGE", names(prep.out$control)[6], sep = ':')
    rownames(varcomp)[
      rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
                                                    grepl('_2', rownames(varcomp)) &
                                                    grepl(names(prep.out$control)[6],
                                                          rownames(varcomp))),])
    ] = paste("IGE", names(prep.out$control)[6], sep = ':')
    # rownames(varcomp)[
    #   rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
    #                                                 grepl('cor', rownames(varcomp)) &
    #                                                 grepl(names(prep.out$control)[6],
    #                                                       rownames(varcomp))),])
    # ] = paste("cor(IGE_DGE)", names(prep.out$control)[6], sep = ':')
    
    if(prep.out$control[,7] > 0){
      rownames(varcomp)[
        rownames(varcomp) %in% rownames(varcomp[which(grepl(paste(names(prep.out$control)[6], 'cor', sep='!'), 
                                                            rownames(varcomp))),])
      ] = paste(paste0("R=autocor(", 
                       paste(names(prep.out$control[7]),
                             1:as.numeric(prep.out$control[7]), sep = '_'), ')'),
                names(prep.out$control)[6], sep = '!')
      
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl(paste(names(prep.out$control)[6], 
                                             '_', sep=''),
                                       rownames(varcomp))),])
      ] = paste0("R=", rownames(varcomp[which(grepl(paste(names(prep.out$control)[6], 
                                                          '_', sep=''),
                                                    rownames(varcomp))),]))
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
        rownames(varcomp) == rownames(varcomp[which(grepl(paste(names(prep.out$control)[6], 'cor', sep='!'), 
                                                          rownames(varcomp))),])
      ] = paste0("R=autocor(", names(prep.out$control)[6],')')
      
      rownames(varcomp)[
        rownames(varcomp) == rownames(varcomp[which(grepl(paste(names(prep.out$control)[6], '_', sep=''), 
                                                          rownames(varcomp))),] )
      ] = paste0("R=", paste(names(prep.out$control)[6],
                             1:as.numeric(prep.out$control[6]), sep = '_'))
      
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('row!',
                                       rownames(varcomp))),])
      ] = 'R=autocor(row)'
      
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('col',
                                       rownames(varcomp))),])
      ] = 'R=autocor(col)'
      
    }
  }else{
    
    if(prep.out$control[,7] > 0){
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
        rownames(varcomp) %in% rownames(varcomp[which(grepl('!R', 
                                                            rownames(varcomp))),])
      ] = 'R'
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('row!',
                                       rownames(varcomp))),])
      ] = 'R=autocor(row):'
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('col',
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
                                                  grepl('cor', rownames(varcomp))),])
  ] = "cor(IGE_DGE)"
  
  output$varcomp = varcomp
  
  
  # BLUPs --------------------
  blup = as.data.frame(summary(model, coef = TRUE)$coef.random)
  
  if(prep.out$control[,6] > 0){
    ## Main effects --------------
    DGE = blup[which(grepl(names(prep.out$control)[2], rownames(blup)) & 
                       !grepl(names(prep.out$control)[6], rownames(blup))), -3]
    DGE[,names(prep.out$control)[2]] = gsub(paste0(names(prep.out$control)[2],'_'), 
                                               '', rownames(DGE), fixed = T)
    rownames(DGE) = NULL
    DGE$rel.DGE = 1-(DGE$std.error^2/varcomp['DGE','component'])
    DGE = DGE[,c(3,1,2,4)]; colnames(DGE)[c(2,3)] = c('DGE', 'se.DGE')
    
    IGE = blup[which(grepl('grp', rownames(blup)) & 
                       !grepl(names(prep.out$control)[6], rownames(blup))), -3]
    IGE[,names(prep.out$control)[2]] = regmatches(rownames(IGE), 
                                                  m = regexpr(paste(DGE$clone, collapse = "|"),
                                                              rownames(IGE)))
    rownames(IGE) = NULL
    IGE$rel.IGE = 1-(IGE$std.error^2/varcomp['IGE','component'])
    IGE = IGE[,c(3,1,2,4)]; colnames(IGE)[c(2,3)] = c('IGE', 'se.IGE')
    
    IGE$class = ifelse(IGE$IGE > mean(IGE$IGE) + stats::sd(IGE$IGE), 'Sensitive', 
                       ifelse(IGE$IGE < mean(IGE$IGE) - stats::sd(IGE$IGE), 'Aggressive', 
                              'Homeostatic'))

    main = merge(DGE, IGE, by = names(prep.out$control)[2])
    
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
    DGE.int = blup[which(grepl(names(prep.out$control)[2], rownames(blup)) & 
                           grepl(names(prep.out$control)[6], rownames(blup))), -3]
    DGE.int[,names(prep.out$control)[2]] = gsub(paste0(names(prep.out$control)[2],'_'),'', 
                                                gsub(':.*', '', rownames(DGE.int)), fixed = T)
    DGE.int[,names(prep.out$control)[6]] = gsub(paste0(names(prep.out$control)[6], '_'),'', 
                                                gsub('.*:', '', rownames(DGE.int)), fixed = T)
    rownames(DGE.int) = NULL
    DGE.int$rel.DGE = 1-(DGE.int$std.error^2/varcomp['DGE','component'])
    DGE.int = DGE.int[,c(3,4,1,2,5)]; colnames(DGE.int)[c(3,4)] = c('DGE', 'se.DGE')
    DGE.int = merge(DGE.int, DGE[,c(names(prep.out$control)[2], 'DGE')],
                    by = names(prep.out$control)[2])
    DGE.int$DGE = DGE.int$DGE.x + DGE.int$DGE.y
    DGE.int = DGE.int[,c(1, 2, 7, 4, 5)]
    
    IGE.int = blup[which(grepl('grp', rownames(blup)) & 
                           grepl(names(prep.out$control)[6], rownames(blup))), -3]
    IGE.int[,names(prep.out$control)[2]] = regmatches(rownames(IGE.int), 
                                                      m = regexpr(paste(DGE$clone, collapse = "|"),
                                                                  rownames(IGE.int)))
    IGE.int[,names(prep.out$control)[6]] = gsub(paste0(names(prep.out$control)[6], '_'),'', 
                                                   gsub('.*:', '', rownames(IGE.int)), fixed = T)
    rownames(IGE.int) = NULL
    IGE.int$rel.IGE = 1-(IGE.int$std.error^2/varcomp['IGE','component'])
    IGE.int = IGE.int[,c(3,4,1,2,5)]; colnames(IGE.int)[c(3,4)] = c('IGE', 'se.IGE')
    IGE.int = merge(IGE.int, IGE[,c(names(prep.out$control)[2], 'IGE')],
                    by = names(prep.out$control)[2])
    IGE.int$IGE = IGE.int$IGE.x + IGE.int$IGE.y
    IGE.int = IGE.int[,c(1, 2, 7, 4, 5)]
    IGE.int$class = do.call(c, lapply(split(IGE.int, IGE.int[,2]), function(x){
      ifelse(x$IGE > mean(x$IGE) + stats::sd(x$IGE), 'Sensitive', 
             ifelse(x$IGE < mean(x$IGE) - stats::sd(x$IGE), 'Aggressive', 
                    'Homeostatic'))
    }))
    
    within = merge(DGE.int, IGE.int, 
                   by = c(names(prep.out$control)[2],
                          names(prep.out$control)[6]))

    if(isFALSE(weight.tgv)){
      within =  merge(within, do.call(rbind, lapply(split(within, within[,2]), function(x){
        data.frame(
          gen = x[,1], 
          age = x[,2],
          TGV = x$DGE + prep.out[[grep(unique(x[,2]),names(prep.out))]]$CIF *
            x$IGE
        )
      })), by.x = c(names(prep.out$control)[2], names(prep.out$control)[6]), 
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
      })), by.x = c(names(prep.out$control)[2], names(prep.out$control)[6]), 
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
      })), by.x = c(names(prep.out$control)[2], names(prep.out$control)[6]), 
      by.y = c('gen','age'))
    }
    
    output$blups = list(main = main, within = within)
    
    ## Other random effects ----------------
    if(any(!grepl("grp\\(g1\\)", attr(model$formulae$random, 'term.labels')))) {
      coef.random = blup[which(!grepl('grp', rownames(blup)) & 
                                 !grepl(names(prep.out$control)[2], rownames(blup))), -3]
      output$blups$coef.random = coef.random
    }
  }else{
    
    ## Main effects --------------
    DGE = blup[which(grepl(names(prep.out$control)[2], rownames(blup)) & 
                       !grepl(names(prep.out$control)[6], rownames(blup))), -3]
    DGE[,names(prep.out$control)[2]] = gsub(paste0(names(prep.out$control)[2],'_'), 
                                               '', rownames(DGE), fixed = T)
    rownames(DGE) = NULL
    DGE$rel.DGE = 1-(DGE$std.error^2/varcomp['DGE','component'])
    DGE = DGE[,c(3,1,2,4)]; colnames(DGE)[c(2,3)] = c('DGE', 'se.DGE')
    
    IGE = blup[which(grepl('grp', rownames(blup)) & 
                       !grepl(names(prep.out$control)[6], rownames(blup))), -3]
    IGE[,names(prep.out$control)[2]] = regmatches(rownames(IGE), 
                                                  m = regexpr(paste(DGE$clone, collapse = "|"), rownames(IGE)))
    rownames(IGE) = NULL
    IGE$rel.IGE = 1-(IGE$std.error^2/varcomp['IGE','component'])
    IGE = IGE[,c(3,1,2,4)]; colnames(IGE)[c(2,3)] = c('IGE', 'se.IGE')
    
    IGE$class = ifelse(IGE$IGE > mean(IGE$IGE) + stats::sd(IGE$IGE), 'Sensitive', 
                       ifelse(IGE$IGE < mean(IGE$IGE) - stats::sd(IGE$IGE), 'Aggressive', 
                              'Homeostatic'))
    
    main = merge(DGE, IGE, by = names(prep.out$control)[2])    
    
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
                                 !grepl(names(prep.out$control)[2], rownames(blup))), -3]
      output$blups = list(main = main, coef.random = coef.random)
    }else{
      output$blups = main
    }
  }
  
  # Plots -----------------------------
  ## Main effects
  output$plots$main = list()
  
  ### Density
  den.ige.main = ggplot(data = main) + 
    geom_density(aes(x = .data$IGE, after_stat(density), fill = 'Homeostatic'), 
                 alpha = .8)
  d = ggplot2::ggplot_build(den.ige.main)$data[[1]]
  den.ige.main = den.ige.main + geom_area(data = subset(d, x > mean(IGE$IGE) + stats::sd(IGE$IGE)), 
                aes(x = .data$x, y = .data$y, fill = 'Sensitive'), alpha = 0.8) + 
    geom_area(data = subset(d, x < mean(IGE$IGE) - stats::sd(IGE$IGE)), 
              aes(x = .data$x, y = .data$y, fill = 'Aggressive'), alpha = .8) +
    geom_density(data = main, aes(x = .data$IGE, after_stat(density)), 
                 color = 'black', linewidth = 1.2, show.legend = F) +
    scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
    labs(x = "IGE", y = 'Density', fill = 'Class') + theme_minimal() +
    theme(legend.position = 'top')
  
  output$plots$main$IGE.density = den.ige.main
  
  ### DGE vs IGE
  plott = ggplot(data = main, aes(x = .data$IGE, y = .data$DGE)) + 
    geom_point(aes(fill = .data$class), colour = 'black', pch = 21, size = 2.2) +
    scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
    labs(fill = 'Class') + theme_minimal() +
    theme(legend.position = 'top') 
  
  output$plots$main$DGEvsIGE = plott
  
  ### DGE and IGE
  temp = data.frame(
    gen = rep(main[, 1], times = 2), 
    blup = c(main$DGE, main$IGE),
    se = c(main$se.DGE, main$se.IGE),
    rel = c(main$rel.DGE, main$rel.IGE),
    comp = rep(c('DGE', 'IGE'), each = prep.out$control[,2])
  )
  
  temp$gen = factor(temp$gen, levels = temp[order(temp$blup, decreasing = TRUE)[
    order(temp$blup, decreasing = TRUE) %in% which(temp$comp == 'DGE')], 1
    ])
  
  plott = ggplot(data = temp, aes(x = .data$gen, y = .data$blup)) + 
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
    labs(x = names(prep.out$control)[2], y = 'BLUP', fill = "Reliability") + 
    guides(fill = guide_colourbar(title.position = 'top', title.hjust = .5, 
                                  barwidth = 9, frame.colour = 'black'))
    
  output$plots$main$DGE.IGE = plott
  
  ### TGV 
  temp = data.frame(gen = main[,1], tgv = main$TGV)
 plott = ggplot(data = temp) + 
    geom_segment(aes(x = stats::reorder(.data$gen, .data$tgv),
                     xend = .data$gen, y = 0, yend = .data$tgv),
                 linewidth = 1) + 
    geom_point(aes(x = stats::reorder(.data$gen, .data$tgv), y = .data$tgv),
               size = 2, color = 'black', alpha = 0.7, shape = 21, stroke = .4, 
               fill = 'darkred') + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, size = 8)) + 
    labs(x = names(prep.out$control)[2], y = 'TGV')
 
 output$plots$main$TGV = plott
  
  if(prep.out$control[,6] > 0) ## Per year ----------------
  {
    output$plots$within = list()
    
    temp = split(within, within[,2])
    
    ## Density
    den.ige.int = do.call(
      ggpubr::ggarrange, 
      c(lapply(temp, function(x){
        aa = ggplot(data = x) + 
          geom_density(aes(x = IGE, after_stat(density), fill = 'Homeostatic'), 
                       alpha = .8)
        d = ggplot2::ggplot_build(aa)$data[[1]]
        aa = aa + geom_area(data = subset(d, x > mean(IGE$IGE) + stats::sd(IGE$IGE)), 
                            aes(x = .data$x, y = .data$y, fill = 'Sensitive'), alpha = 0.8) + 
          geom_area(data = subset(d, x < mean(IGE$IGE) - stats::sd(IGE$IGE)), 
                    aes(x = .data$x, y = .data$y, fill = 'Aggressive'), alpha = .8) +
          geom_density(data = x, aes(x = .data$IGE, after_stat(density)), 
                       color = 'black', linewidth = 1.2, show.legend = F) +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
          labs(x = "IGE", y = 'Density', fill = 'Class', 
               subtitle = paste(names(prep.out$control)[6], unique(x[,2]))) + 
          theme_minimal() +
          theme(legend.position = 'top')
      }), common.legend = TRUE)
    )
    
    output$plots$within$IGE.density = den.ige.int
    
    ### DGE vs IGE
    plott = do.call(
      ggpubr::ggarrange, 
      c(lapply(temp, function(x){
        aa =  ggplot(data = x, aes(x = .data$IGE, y = .data$DGE)) + 
          geom_point(aes(fill = .data$class), colour = 'black', pch = 21, size = 2.2) +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
          labs(fill = 'Class', subtitle = paste(names(prep.out$control)[6],
                                                unique(x[,2]))) +
          theme_minimal() +
          theme(legend.position = 'top') 
      }), common.legend = TRUE)
    )
    
    output$plots$within$DGEvsIGE = plott
    
    ### DGE and IGE
    plott = do.call(
      ggpubr::ggarrange, 
      c(lapply(temp, function(x){
        
        temp2 = data.frame(
          gen = rep(x[, 1], times = 2), 
          blup = c(x$DGE, x$IGE),
          se = c(x$se.DGE, x$se.IGE),
          rel = c(x$rel.DGE, x$rel.IGE),
          comp = rep(c('DGE', 'IGE'), each = prep.out$control[,2])
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
          theme(axis.text.x = element_text(angle = 90, size = 8), 
                legend.position = 'top', panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey')) +
          scale_fill_viridis_c(option = 'turbo') + 
          labs(x = names(prep.out$control)[2], y = 'BLUP', fill = "Reliability") + 
          guides(fill = guide_colourbar(title.position = 'top', title.hjust = .5, 
                                        barwidth = 9, frame.colour = 'black'))
        
      }), common.legend = TRUE)
    )
    
    output$plots$within$DGE.IGE = plott
    
    ### TGV 
    plott = do.call(
      ggpubr::ggarrange, 
      c(lapply(temp, function(x){
        
        aa = data.frame(gen = x[,1], tgv = x$TGV)
        
        bb = ggplot(data = aa) + 
          geom_segment(aes(x = stats::reorder(.data$gen, .data$tgv),
                           xend = .data$gen, y = 0, yend = .data$tgv),
                       linewidth = 1) + 
          geom_point(aes(x = stats::reorder(.data$gen, .data$tgv), y = .data$tgv),
                     size = 2, color = 'black', alpha = 0.7, shape = 21, stroke = .4, 
                     fill = 'darkred') + 
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, size = 8)) + 
          labs(x = names(prep.out$control)[2], y = 'TGV',
               subtitle = paste(names(prep.out$control)[6],
                                unique(x[,2])))
      }))
    )
    
    output$plots$within$TGV = plott

    ## Number of neighbours
    temp = split(prep.out$data, prep.out$data[,names(prep.out$control)[6]])
    
    nneigh = do.call(
      ggpubr::ggarrange,
      c(lapply(temp, function(x){
        Z = x[, 1:prep.out$control[,2]]
        Z = ifelse(crossprod(ifelse(Z == 0, 0, 1)) == 0, 0, 1)
        diag(Z) = 0
        Z = Z[order(rownames(Z)), order(colnames(Z))]
        
        neigh = do.call(rbind, lapply(
          apply(Z, 1, function(y){
            data.frame(
              neigh = names(which(y == 1)),
              class = within[
                which(within[,2] == unique(x[, names(prep.out$control)[6]])),
                'class'
              ][
                match(names(which(y == 1)),
                      within[which(within[,2] == unique(x[, names(prep.out$control)[6]])),
                             c(names(prep.out$control)[2])])
              ]
            )
          }), function(w){
            data.frame(table(w$class))
          }
        ))
        
        neigh$gen = substr(rownames(neigh), start = 1, stop = nchar(unique(rownames(Z))))
        
        ggplot(data = neigh, aes(x = .data$gen, y = .data$Freq, fill = .data$Var1)) +
          geom_bar(stat = "identity", color = 'black') +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, size = 8),
                axis.title = element_text(face = 'bold'),
                legend.position = 'bottom') +
          labs(x = names(prep.out$control)[2], 
               y = "No of different genotypes as neighbours", fill = 'Class', 
               subtitle = paste(names(prep.out$control)[6], 
                                unique(x[,names(prep.out$control)[6]]))) +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6'))
      }), common.legend = TRUE)
    )
    
    output$plots$within$n.neigh = nneigh
    
    ### Grids 
    if(prep.out$control[,7] > 0){
      dat = as.data.frame(model$mf)[,names(prep.out$control)[c(2,4,5,7,6,1)]]
      dat$e = c(model$residuals)
      colnames(dat) = c('gen', 'row', 'col', 'area', 'age', 'y', 'e')
      temp = merge(dat, within[,c(1,2,3,6,9)], by.x = c('gen', 'age'), 
                   by.y = c(names(prep.out$control)[2],
                            names(prep.out$control)[6]))
      
      facet.label.col =  paste(names(prep.out$control)[7],
                           1:as.numeric(prep.out$control[7]))
      names(facet.label.col) = unique(temp$area)
      facet.label.row =  paste(names(prep.out$control)[6],
                               1:as.numeric(prep.out$control[6]))
      names(facet.label.row) = unique(temp$age)
      
      output$plots$within$grid.pheno = ggplot(data = temp, 
                                            aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                fill = .data$y)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), cols = vars(.data$area), 
                   scales = 'free_x', 
                   labeller = labeller(.cols = facet.label.col,
                                       .rows = facet.label.row)) +
        scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'Phenotype') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
      output$plots$within$grid.DGE = ggplot(data = temp, 
                                             aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                 fill = .data$DGE)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), cols = vars(.data$area), 
                   scales = 'free_x', 
                   labeller = labeller(.cols = facet.label.col,
                                       .rows = facet.label.row)) +
        scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'DGE') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
      output$plots$within$grid.IGE = ggplot(data = temp, 
                                            aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                fill = .data$IGE)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), cols = vars(.data$area), 
                   scales = 'free_x', 
                   labeller = labeller(.cols = facet.label.col,
                                       .rows = facet.label.row)) +
        scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'IGE') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
      output$plots$within$grid.IGE.class = ggplot(data = temp, 
                                                  aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                      fill = .data$class)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), cols = vars(.data$area), 
                   scales = 'free_x', 
                   labeller = labeller(.cols = facet.label.col,
                                       .rows = facet.label.row)) +
        scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'IGE: Classes') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'bottom', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
    }else{
      dat = as.data.frame(model$mf)[,names(prep.out$control)[c(2,4,5,6,1)]]
      dat$e = c(model$residuals)
      colnames(dat) = c('gen', 'row', 'col', 'age', 'y', 'e')
      temp = merge(dat, within[,c(1,2,3,6,9)], by.x = c('gen', 'age'), 
                   by.y = c(names(prep.out$control)[2],
                            names(prep.out$control)[6]))
      
      facet.label.row =  paste(names(prep.out$control)[6],
                               1:as.numeric(prep.out$control[6]))
      names(facet.label.row) = unique(temp$age)
      
      output$plots$within$grid.pheno = ggplot(data = temp, 
                                              aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                  fill = .data$y)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), scales = 'free_x', 
                   labeller = labeller(.rows = facet.label.row)) +
        scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'Phenotype') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
      output$plots$within$grid.DGE = ggplot(data = temp, 
                                            aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                fill = .data$DGE)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), scales = 'free_x', 
                   labeller = labeller(.rows = facet.label.row)) +
        scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'DGE') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
      output$plots$within$grid.IGE = ggplot(data = temp, 
                                            aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                fill = .data$IGE)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), scales = 'free_x', 
                   labeller = labeller(.rows = facet.label.row)) +
        scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'IGE') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
      output$plots$within$grid.IGE.class = ggplot(data = temp, 
                                                  aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                      fill = .data$class)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), scales = 'free_x', 
                   labeller = labeller(.rows = facet.label.row)) +
        scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'IGE: Classes') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'bottom', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
    }
    
  }else{
    ## Number of neighbours
    Z = prep.out$Z
    Z = ifelse(crossprod(ifelse(Z == 0, 0, 1)) == 0, 0, 1)
    diag(Z) = 0
    Z = Z[order(rownames(Z)), order(colnames(Z))]
    
    neigh = do.call(rbind, lapply(
      apply(Z, 1, function(y){
        data.frame(
          neigh = names(which(y == 1)),
          class = main[match(names(which(y == 1)), main[,1]), 8]
        )
      }), function(w){
        data.frame(table(w$class))
      }
    ))
    
    neigh$gen = substr(rownames(neigh), start = 1, stop = nchar(unique(rownames(Z))))
    
    nneigh = ggplot(data = neigh, aes(x = .data$gen, y = .data$Freq, 
                                      fill = .data$Var1)) +
      geom_bar(stat = "identity", color = 'black') +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, size = 8),
            axis.title = element_text(face = 'bold'),
            legend.position = 'bottom') +
      labs(x = names(prep.out$control)[2], 
           y = "Number of neighbours", fill = 'Class') +
      scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6'));nneigh
    
    output$plots$main$n.neigh = nneigh
    
    
    ### Grids 
    if(prep.out$control[,7] > 0){
      
      dat = as.data.frame(model$mf)[,names(prep.out$control)[c(2,4,5,7,1)]]
      dat$e = c(model$residuals)
      colnames(dat) = c('gen', 'row', 'col', 'area', 'y', 'e')
      temp = merge(dat, main[,c(1,2,5,8)], by.x = 'gen', by.y = names(prep.out$control)[2])
      
      facet.label =  paste(names(prep.out$control)[7],
                           1:as.numeric(prep.out$control[7]))
      names(facet.label) = unique(temp$area)
      
      plott = ggpubr::ggarrange(
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$y)) + 
          geom_tile(color = 'black') + 
          facet_wrap(.~.data$area, scales = 'free_x', 
                     labeller = labeller(.cols = facet.label)) +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = '',
               title = 'Phenotype') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey')),
        # ggplot(data = temp, 
        #        aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
        #            fill = .data$e)) + 
        #   geom_tile(color = 'black') + 
        #   facet_wrap(.~.data$area, scales = 'free_x', 
        #              labeller = labeller(.cols = facet.label)) +
        #   scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        #   labs(x = "Row", y = 'Column', fill = '',
        #        title = 'Residue') + 
        #   theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
        #         panel.background = element_blank(), 
        #         panel.grid = element_line(colour = 'lightgrey')), 
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$DGE)) + 
          geom_tile(color = 'black') + 
          facet_wrap(.~.data$area, scales = 'free_x', 
                     labeller = labeller(.cols = facet.label)) +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = '',
               title = 'DGE') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey')),
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$class)) + 
          geom_tile(color = 'black') + 
          facet_wrap(.~.data$area, scales = 'free_x', 
                     labeller = labeller(.cols = facet.label)) +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
          labs(x = "Row", y = 'Column', fill = '',
               title = 'IGE: Classes') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'top', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey')),
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$IGE)) + 
          geom_tile(color = 'black') + 
          facet_wrap(.~.data$area, scales = 'free_x', 
                     labeller = labeller(.cols = facet.label)) +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = '',
               title = 'IGE') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))
      )
      
      output$plots$main$grid = plott
      
    }else{
      temp = merge(
        data.frame(prep.out$neigh_check[,1:4], e = model$residuals),
        main[,c(1,2,5,8)], by.x = 'gen', by.y = names(prep.out$control)[2]
      )
      
      plott = ggpubr::ggarrange(
        ggplot(data = temp, 
               aes(x = .data$row, y = .data$col, fill = .data$y_focal)) + 
          geom_tile(color = 'black') +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = '',
               title = 'Phenotype') + 
          theme_minimal() + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right'),
        # ggplot(data = temp, 
        #        aes(x = .data$row, y = .data$col, fill = .data$e)) + 
        #   geom_tile(color = 'black') + 
        #   scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        #   labs(x = "Row", y = 'Column', fill = '', title = 'Residual') + 
        #   theme_minimal()+ 
        #   theme(plot.title = element_text(hjust = .5)), 
        ggplot(data = temp, 
               aes(x = .data$row, y = .data$col, fill = .data$DGE)) + 
          geom_tile(color = 'black') + 
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = '', title = 'DGE') + 
          theme_minimal()+ 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right'),
        ggplot(data = temp, 
               aes(x = .data$row, y = .data$col, fill = .data$class)) + 
          geom_tile(color = 'black')+
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
          labs(x = "Row", y = 'Column', fill = '', title = 'IGE: Classes') + 
          theme_minimal()+ 
          theme(plot.title = element_text(hjust = .5), legend.position = 'top'),
        ggplot(data = temp, 
               aes(x = .data$row, y = .data$col, fill = .data$IGE)) + 
          geom_tile(color = 'black') + 
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = '', title = 'IGE') + 
          theme_minimal()+ 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right')
      )
      
      output$plots$main$grid = plott
    }
  }
  remove(prep.out, envir = .GlobalEnv)
  
  class(output) = 'comp.resp'
  
  return(output)
}









