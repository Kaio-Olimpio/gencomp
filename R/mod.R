##' @title Fit a genetic-spatial competition model
##' 
##' @description
##' This function fits genetic-spatial competition models using [asreml::asreml()]
##' 
##' @param prep.out A `comprep` object.
##' @inheritParams asreml::asreml 
##' @param random The function is set to fit the genetic competition models in the 
##' `random` term. To consider further random effects, `random` receives a formula
##'  object specifying them. Otherwise, `random = ~1` (default). This argument has 
##'  the same general characteristics as fixed, but there can be no left side to 
##'  the ~ operator. Variance structures imposed on random terms are specified using 
##'  special model functions. See [asreml::asreml()] for details.
##' @param cor A logical value. If `TRUE` (default), fits a model considering the correlation 
##' between direct and indirect genetic effects.
##' @param lrtest A logical value. If `TRUE`, performs a likelihood ratio test to verify
##' the significance of the direct and indirect genetic effects. Defaults to `FALSE`.
##' @inheritDotParams asreml::asreml -fixed -random -data -residual
##' 
##' @return An object of class `asreml` containing the results of the fitted linear model.
##' Instances of generic methods such as `plot()`, `predict()` and `summary()` return
##' various derived results of the fit. The method `resid()`, `coef()` and `fitted()` 
##' extract some of its components. See [asreml::asreml.object()] for details of the 
##' components of the returned list.
##'
##' @details
##' A general genetic-spatial linear mixed model can be represented by: 
##' 
##' \deqn{\mathbf{y} = \mathbf{X} \boldsymbol{\beta} + \mathbf{Z}_g \mathbf{g} + \mathbf{Z}_c \mathbf{c} + \mathbf{Z}_p \mathbf{p}  + \boldsymbol{\xi}}
##' 
##' where \eqn{\mathbf{y}} is the vector of phenotypic records, \eqn{\boldsymbol{\beta}}
##' is the vector of fixed effects, \eqn{\mathbf{g}} is the vector of direct genetic 
##' effects (DGE), \eqn{\mathbf{c}} is the vector of indirect genetic effects (IGE),
##' \eqn{\mathbf{p}} is the vector of other random effects, 
##' and \eqn{\boldsymbol{\xi}} is the vector of spatially correlated errors. 
##' \eqn{\mathbf{X}} is the incidence matrix of the fixed effects, \eqn{\mathbf{Z}_g} 
##' is the DGE incidence matrix, \eqn{\mathbf{Z}_c} is the IGE incidence matrix (built 
##' using [GenComp::prep]), and \eqn{\mathbf{Z}_p} is the design matrix of 
##' other random effects. The dimensions of \eqn{\mathbf{Z}_c} are the same 
##' as \eqn{\mathbf{Z}_g}. The spatially correlated errors are distributed as 
##' \eqn{\boldsymbol{\xi} \sim N\{\mathbf{0}, \sigma^2_\xi[\mathbf{AR1}(\rho_C) \otimes \mathbf{AR1}(\rho_R)]\}}, 
##' where \eqn{\sigma^2_\xi} is the spatially correlated residual variance, 
##' \eqn{\mathbf{AR1}(\rho_C)} and \eqn{\mathbf{AR1}(\rho_R)} are the first-order autoregressive
##'  correlation matrices in the column and row directions, and \eqn{\otimes} is the
##' Kronecker product. If `cor = TRUE`, the function will fit a model in which
##'  \eqn{\mathbf{g}} and \eqn{\mathbf{c}} are correlated outcomes of the genotypic effects
##'  decomposition. They both follow a Gaussian distribution, with mean centred in zero, and
##'  covariange given by:
##'  
##'  \deqn{\mathbf{\Sigma_g}=  \begin{bmatrix}\sigma_{\text{g}}^2 & \sigma_{\text{gc}}\\\sigma_{\text{gc}} & \sigma_{\text{c}}^2\\\end{bmatrix}\otimes {{\mathbf I_V}}}
##'  
##'  where \eqn{\sigma_{\text{g}}^2} is the DGE variance, \eqn{\sigma_{\text{c}}^2} is the IGE variance, 
##'  and \eqn{\sigma_{\text{gc}}} is the covariance between DGE and IGE.
##'  
##'  If the dataset has multiple ages, the model described aboved is expanded with two more effects:
##'  the DGE \eqn{\times} age interaction, and the IGE \eqn{\times} age interaction.
##'  
##'  The likelihood ratio test is performed using a model without the correlation between DGE and IGE.
##'  
##'
##' @seealso  [GenComp::prep], [asreml::asreml.options], [asreml::asreml.object], [asreml::family_dist]
##'
##' @import asreml
##' @importFrom stats as.formula
##' 
##' @export
##' 
##' @examples
##' \donttest{
##'  library(GenComp)
##'  comp_mat = prep(data = euca, gen = 'clone', repl = 'block', area = 'area', 
##'                       ind = 'tree', age = 'age', row = 'row', col = 'col', 
##'                       dist.col = 3, dist.row = 2, trait = 'mai', method = 'SK',
##'                       n.dec = 3, verbose = TRUE)
##'  
##'  model = asr(prep.out = comp_mat, 
##'              fixed = mai~ age, 
##'              random = ~ block:age, 
##'              cor = TRUE, maxit = 50,
##'              lrtest = FALSE)
##'  }

asr = function(prep.out, fixed, random = ~1, cor = TRUE, lrtest = FALSE,...) {
  
  requireNamespace('asreml')
  
  control = attr(prep.out, 'control')
  
  prep.out <<- prep.out
  control <<- control
  
  input = prep.out$data
  
  if(control[,7] > 0) ### Multi-areas ------------
    {
      if(control[,6] > 0) ### Multi-ages --------
      {
        input = input[order(input[,names(control)[7]], input[,names(control)[6]],
                            input[,names(control)[4]], input[,names(control)[5]]),]
        
        if(cor){
          scm = asreml::asreml(fixed = fixed,
                       random = stats::as.formula(
                         paste(paste(as.character(random), collapse = ""), 
                               paste0('+ str(~', colnames(control)[2],
                                      '+ grp(g1), ~corh(2):id(',
                                      colnames(control)[2],')) + str(~',
                                      colnames(control)[2],':',
                                      colnames(control)[6],'+grp(g1):',
                                      colnames(control)[6], ', ~diag(2):id(',
                                      colnames(control)[6], '):id(',
                                      colnames(control)[2], '))'))
                       ),
                       group = list(g1 = 1:control[,2]),
                       residual = stats::as.formula(paste0(
                         '~dsum(~ar1h(', colnames(control)[6],'):ar1(',
                         colnames(control)[4],'):ar1(',
                         colnames(control)[5],')|',
                         colnames(control)[7],')')
                       ),
                       data = input, ...)

        }else{
          scm = asreml::asreml(fixed = fixed,
                       random = stats::as.formula(
                         paste(paste(as.character(random), collapse = ""), 
                               paste0('+', colnames(control)[2],
                                      '+ grp(g1) +',
                                      colnames(control)[2],':',
                                      colnames(control)[6],'+grp(g1):',
                                      colnames(control)[6]))
                       ),
                       group = list(g1 = 1:control[,2]),
                       residual = stats::as.formula(paste0(
                         '~dsum(~ar1h(', colnames(control)[6],'):ar1(',
                         colnames(control)[4],'):ar1(',
                         colnames(control)[5],')|',
                         colnames(control)[7],')')
                       ),
                       data = input, ...)

        }
      }else ### Single age ---------------
      {
        input = input[order(input[,names(control)[7]],
                            input[,names(control)[4]], 
                            input[,names(control)[5]]),]
        if(cor){
          scm = asreml::asreml(fixed = fixed,
                       random = stats::as.formula(
                         paste(paste(as.character(random), collapse = ""), 
                               paste0('+ str(~', colnames(control)[2],
                                      '+ grp(g1), ~corh(2):id(',
                                      colnames(control)[2],'))'))
                       ),
                       group = list(g1 = 1:control[,2]),
                       residual = stats::as.formula(paste0(
                         '~dsum(~ar1(', colnames(control)[4],'):ar1(',
                         colnames(control)[5],')|',
                         colnames(control)[7],')')
                       ),
                       data = input, ...)
        }else{
          scm = asreml::asreml(fixed = fixed,
                       random = stats::as.formula(
                         paste(paste(as.character(random), collapse = ""), 
                               paste0('+', colnames(control)[2],
                                      '+ grp(g1)'))
                       ),
                       group = list(g1 = 1:control[,2]),
                       residual = stats::as.formula(paste0(
                         '~dsum(~ar1(', colnames(control)[4],'):ar1(',
                         colnames(control)[5],')|',
                         colnames(control)[7],')')
                       ),
                       data = input, ...)
        }
      } 
    }else ### Single area -----------
      {
        if(control[,6] > 0) ### Multi-ages --------
        {
          input = input[order(input[,names(control)[6]],
                              input[,names(control)[4]], 
                              input[,names(control)[5]]),]
          if(cor){
            scm = asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+ str(~', colnames(control)[2],
                                        '+ grp(g1), ~corh(2):id(',
                                        colnames(control)[2],')) + str(~',
                                        colnames(control)[2],':',
                                        colnames(control)[6],'+grp(g1):',
                                        colnames(control)[6], ', ~diag(2):id(',
                                        colnames(control)[6], '):id(',
                                        colnames(control)[2], '))'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1h(', colnames(control)[6],'):ar1(',
                           colnames(control)[4],'):ar1(',
                           colnames(control)[5],')')
                         ),
                         data = input, ...)
            
          }else{
            scm = asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2],
                                        '+ grp(g1) +',
                                        colnames(control)[2],':',
                                        colnames(control)[6],'+grp(g1):',
                                        colnames(control)[6]))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1h(', colnames(control)[6],'):ar1(',
                           colnames(control)[4],'):ar1(',
                           colnames(control)[5],')')
                         ),
                         data = input, ...)
          }
        }else ### Single age ---------------
        {
          input = input[order(input[,names(control)[4]], 
                              input[,names(control)[5]]),]
          if(cor){
            scm = asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+ str(~', colnames(control)[2],
                                        '+ grp(g1), ~corh(2):id(',
                                        colnames(control)[2],'))'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1(', colnames(control)[4],'):ar1(',
                           colnames(control)[5],')')
                         ),
                         data = input, ...)
          }else{
            scm = asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2],
                                        '+ grp(g1)'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1(', colnames(control)[4],'):ar1(',
                           colnames(control)[5],')')
                         ),
                         data = input, ...)
          }
        }
      }
  
  # LRT ---------------------
  if(lrtest){
    message("====> Starting likelihood ratio tests")
    if(control[,7] > 0) ### Multi-areas ------------
    {
      if(control[,6] > 0) ### Multi-ages --------
      {
        input = input[order(input[,names(control)[7]], input[,names(control)[6]],
                            input[,names(control)[4]], input[,names(control)[5]]),]
        
        bench = suppressWarnings({asreml::asreml(fixed = fixed,
                                                 random = stats::as.formula(
                                                   paste(paste(as.character(random), collapse = ""), 
                                                         paste0('+', colnames(control)[2],
                                                                '+ grp(g1) +',
                                                                colnames(control)[2],':',
                                                                colnames(control)[6],'+grp(g1):',
                                                                colnames(control)[6]))
                                                 ),
                                                 group = list(g1 = 1:control[,2]),
                                                 residual = stats::as.formula(paste0(
                                                   '~dsum(~ar1h(', colnames(control)[6],'):ar1(',
                                                   colnames(control)[4],'):ar1(',
                                                   colnames(control)[5],')|',
                                                   colnames(control)[7],')')
                                                 ),
                                                 data = input,...)})
        
        lr.dint = suppressWarnings({lrt(bench, 
                                        asreml::asreml(fixed = fixed,
                                                       random = stats::as.formula(
                                                         paste(paste(as.character(random), collapse = ""), 
                                                               paste0('+', colnames(control)[2],
                                                                      '+ grp(g1)','+grp(g1):',
                                                                      colnames(control)[6]))
                                                       ),
                                                       group = list(g1 = 1:control[,2]),
                                                       residual = stats::as.formula(paste0(
                                                         '~dsum(~ar1h(', colnames(control)[6],'):ar1(',
                                                         colnames(control)[4],'):ar1(',
                                                         colnames(control)[5],')|',
                                                         colnames(control)[7],')')
                                                       ),
                                                       data = input,...))})
        
        lr.iint = suppressWarnings({
          lrt(bench, 
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+', colnames(control)[2],
                                            '+ grp(g1) +',
                                            colnames(control)[2],':',
                                            colnames(control)[6]))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~dsum(~ar1h(', colnames(control)[6],'):ar1(',
                               colnames(control)[4],'):ar1(',
                               colnames(control)[5],')|',
                               colnames(control)[7],')')
                             ),
                             data = input, ...))
        }) 
        
        lr.d = suppressWarnings({
          lrt(bench, 
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+ grp(g1)','+grp(g1):',
                                            colnames(control)[6]))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~dsum(~ar1h(', colnames(control)[6],'):ar1(',
                               colnames(control)[4],'):ar1(',
                               colnames(control)[5],')|',
                               colnames(control)[7],')')
                             ),
                             data = input,...))
        }) 
        
        lr.i = suppressWarnings({
          lrt(bench, 
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+', colnames(control)[2],
                                            '+',
                                            colnames(control)[2],':',
                                            colnames(control)[6]))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~dsum(~ar1h(', colnames(control)[6],'):ar1(',
                               colnames(control)[4],'):ar1(',
                               colnames(control)[5],')|',
                               colnames(control)[7],')')
                             ),
                             data = input,...))
        }) 
        
        lr.result = cbind(effect = c("DGE", 'IGE', 'DGE:age', 'IGE:age'), 
                          rbind(lr.d[,-1],lr.i[,-1], lr.dint[,-1], lr.iint[,-1]))
        rownames(lr.result) = NULL
        
      }else ### Single age ---------------
      {
        input = input[order(input[,names(control)[7]],
                            input[,names(control)[4]], input[,names(control)[5]]),]
        
        bench = suppressWarnings({
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2],
                                        '+ grp(g1)'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~dsum(~ar1(', colnames(control)[4],'):ar1(',
                           colnames(control)[5],')|',
                           colnames(control)[7],')')
                         ),
                         data = input, ...)
        }) 
        
        lr.d = suppressWarnings({
          lrt(bench,  asreml::asreml(fixed = fixed,
                                     random = stats::as.formula(
                                       paste(paste(as.character(random), collapse = ""), 
                                             paste0('+ grp(g1)'))
                                     ),
                                     group = list(g1 = 1:control[,2]),
                                     residual = stats::as.formula(paste0(
                                       '~dsum(~ar1(', colnames(control)[4],'):ar1(',
                                       colnames(control)[5],')|',
                                       colnames(control)[7],')')
                                     ),
                                     data = input, ...))
        }) 
        
        lr.i = suppressWarnings({
          lrt(bench, asreml::asreml(fixed = fixed,
                                    random = stats::as.formula(
                                      paste(paste(as.character(random), collapse = ""), 
                                            paste0('+', colnames(control)[2]))
                                    ),
                                    group = list(g1 = 1:control[,2]),
                                    residual = stats::as.formula(paste0(
                                      '~dsum(~ar1(', colnames(control)[4],'):ar1(',
                                      colnames(control)[5],')|',
                                      colnames(control)[7],')')
                                    ),
                                    data = input, ...))
        }) 
        lr.result = cbind(effect = c("DGE", 'IGE'), 
                          rbind(lr.d[,-1],lr.i[,-1]))
        rownames(lr.result) = NULL
      } 
    }else ### Single area -----------
    {
      if(control[,6] > 0) ### Multi-ages --------
      {
        input = input[order(input[,names(control)[6]],
                            input[,names(control)[4]], input[,names(control)[5]]),]
        bench = suppressWarnings({
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2],
                                        '+ grp(g1) +',
                                        colnames(control)[2],':',
                                        colnames(control)[6],'+grp(g1):',
                                        colnames(control)[6]))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1h(', colnames(control)[6],'):ar1(',
                           colnames(control)[4],'):ar1(',
                           colnames(control)[5],')')
                         ),
                         data = input, ...)
        }) 
        
        lr.dint = suppressWarnings({
          lrt(bench, 
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+', colnames(control)[2],
                                            '+ grp(g1)','+grp(g1):',
                                            colnames(control)[6]))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~ar1h(', colnames(control)[6],'):ar1(',
                               colnames(control)[4],'):ar1(',
                               colnames(control)[5],')')
                             ),
                             data = input, ...))
        }) 
        
        lr.iint = suppressWarnings({
          lrt(bench, 
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+', colnames(control)[2],
                                            '+ grp(g1) +',
                                            colnames(control)[2],':',
                                            colnames(control)[6]))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~ar1h(', colnames(control)[6],'):ar1(',
                               colnames(control)[4],'):ar1(',
                               colnames(control)[5],')')
                             ),
                             data = input, ...))
        }) 
        
        lr.d = suppressWarnings({
          lrt(bench, 
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+ grp(g1) + grp(g1):',
                                            colnames(control)[6]))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~ar1h(', colnames(control)[6],'):ar1(',
                               colnames(control)[4],'):ar1(',
                               colnames(control)[5],')')
                             ),
                             data = input, ...))
        })
        
        lr.i = suppressWarnings({
          lrt(bench, 
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+', colnames(control)[2],
                                            ' +',
                                            colnames(control)[2],':',
                                            colnames(control)[6]))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~ar1h(', colnames(control)[6],'):ar1(',
                               colnames(control)[4],'):ar1(',
                               colnames(control)[5],')')
                             ),
                             data = input, ...))
        }) 
        
        lr.result = cbind(effect = c("DGE", 'IGE', 'DGE:age', 'IGE:age'), 
                          rbind(lr.d[,-1],lr.i[,-1], lr.dint[,-1], lr.iint[,-1]))
        rownames(lr.result) = NULL
        
        
      }else ### Single age ---------------
      {
        input = input[order(input[,names(control)[4]], input[,names(control)[5]]),]
        
        bench = suppressWarnings({
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2],
                                        '+ grp(g1)'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1(', colnames(control)[4],'):ar1(',
                           colnames(control)[5],')')
                         ),
                         data = input, ...)
        }) 
        
        lr.d = suppressWarnings({
          lrt(bench, 
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+ grp(g1)'))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~ar1(', colnames(control)[4],'):ar1(',
                               colnames(control)[5],')')
                             ),
                             data = input, ...))
        }) 
        
        lr.i = suppressWarnings({
          lrt(bench, 
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+', colnames(control)[2]))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~ar1(', colnames(control)[4],'):ar1(',
                               colnames(control)[5],')')
                             ),
                             data = input, ...))
        }) 
        
        lr.result = cbind(effect = c("DGE", 'IGE'), 
                          rbind(lr.d[,-1],lr.i[,-1]))
        rownames(lr.result) = NULL
        
      }
    }
    scm$lrt = lr.result
  }

  remove(prep.out, envir = .GlobalEnv)
  remove(control, envir = .GlobalEnv)

  
  return(scm)
}


