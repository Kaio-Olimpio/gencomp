# Single-environment spatial competition model ----------------------------

##' @title Fit a genetic competition model
##' 
##' @description
##' This function fits genetic competition models using [asreml::asreml()]
##' 
##' @param prep.out A `comprep` object.
##' @inheritParams asreml::asreml 
##' @param random The function is set to fit the genetic competition models in the 
##' `random` term. To consider further random effects, `random` receives a formula
##'  object specifying them. Otherwise, `random = ~1` (default). This argument has 
##'  the same general characteristics as fixed, but there can be no left side to 
##'  the ~ operator. Variance structures imposed on random terms are specified using 
##'  special model functions. See [asreml::asreml()] for details.
##' @param spatial A logical value. If `TRUE` (default), fits a spatial-genetic competition model 
##' @param cor A logical value. If `TRUE` (default), fits a model considering the correlation 
##' between direct and indirect genetic effects.
##' @param lrtest A logical value. If `TRUE`, performs a likelihood ratio test to verify
##' the significance of the direct and indirect genetic effects. Defaults to `FALSE`.
##' @inheritDotParams asreml::asreml -fixed -random -data -residual -group
##' 
##' @return An object of class `asreml` containing the results of the fitted linear model.
##' Instances of generic methods such as `plot()`, `predict()` and `summary()` return
##' various derived results of the fit. The method `resid()`, `coef()` and `fitted()` 
##' extract some of its components. See [asreml::asreml.object()] for details of the 
##' components of the returned list.
##'
##' @details
##' A general genetic competition linear mixed model can be represented by: 
##' 
##' \deqn{\mathbf{y} = \mathbf{X} \boldsymbol{\beta} + \mathbf{Z}_g \mathbf{g} + \mathbf{Z}_c \mathbf{c} + \mathbf{Z}_p \mathbf{p}  + \boldsymbol{\varepsilon}}
##' 
##' where \eqn{\mathbf{y}} is the vector of phenotypic records, \eqn{\boldsymbol{\beta}}
##' is the vector of fixed effects, \eqn{\mathbf{g}} is the vector of direct genetic 
##' effects (DGE), \eqn{\mathbf{c}} is the vector of indirect genetic effects (IGE),
##' \eqn{\mathbf{p}} is the vector of other random effects, 
##' and \eqn{\boldsymbol{\varepsilon}} is the vector of errors.
##' \eqn{\mathbf{X}} is the incidence matrix of the fixed effects, \eqn{\mathbf{Z}_g} 
##' is the DGE incidence matrix, \eqn{\mathbf{Z}_c} is the IGE incidence matrix (the genetic competition matrix), 
##' and \eqn{\mathbf{Z}_p} is the design matrix of 
##' other random effects. The dimensions of \eqn{\mathbf{Z}_c} are the same 
##' as \eqn{\mathbf{Z}_g}. If `spatial = TRUE`, \eqn{\varepsilon} is a vector of spatially correlated errors, distributed as 
##' \eqn{\boldsymbol{\varepsilon} \sim N\{\mathbf{0}, \sigma^2_\varepsilon[\mathbf{AR1}(\rho_C) \otimes \mathbf{AR1}(\rho_R)]\}}, 
##' where \eqn{\sigma^2_\varepsilon} is the spatially correlated residual variance, 
##' \eqn{\mathbf{AR1}(\rho_C)} and \eqn{\mathbf{AR1}(\rho_R)} are the first-order autoregressive
##'  correlation matrices in the column and row directions, and \eqn{\otimes} is the
##' Kronecker product. If `cor = TRUE`, the function will fit a model in which
##'  \eqn{\mathbf{g}} and \eqn{\mathbf{c}} are correlated outcomes of the genotypic effects'
##'  decomposition. They both follow a Gaussian distribution, with mean centred in zero, and
##'  covariance given by:
##'  
##'  \deqn{\mathbf{\Sigma_g} = \begin{bmatrix}\sigma_{\text{g}}^2 & \sigma_{\text{gc}}\\\sigma_{\text{gc}} & \sigma_{\text{c}}^2\\\end{bmatrix}\otimes {{\mathbf I_V}}}
##'  
##'  where \eqn{\sigma_{\text{g}}^2} is the DGE variance, \eqn{\sigma_{\text{c}}^2} is the IGE variance, 
##'  and \eqn{\sigma_{\text{gc}}} is the covariance between DGE and IGE.
##'  
##'  The likelihood ratio test is performed using a model without the correlation between DGE and IGE.
##'  
##'
##' @seealso  [gencomp::prepfor], [gencomp::prepcrop], [asreml::asreml.options], [asreml::asreml.object], [asreml::family_dist]
##'
##' @import asreml
##' @importFrom stats as.formula
##' 
##' @export
##' 
##' @examples
##' \donttest{
##'  library(gencomp)
##'  dat = euca[which(euca$age == "6y"),]
##'  comp_mat = prepfor(data = dat, gen = 'clone', area = 'area',
##'                     ind = 'tree', age = NULL, row = 'row', col = 'col', 
##'                     dist.col = 3, dist.row = 2, trait = 'MAI', method = 'SK',
##'                     n.dec = 3, verbose = FALSE, effs = c("block"))
##'  
##'  model = asr(prep.out = comp_mat, 
##'              fixed = MAI ~ 1, 
##'              random = ~ block,
##'              cor = TRUE,
##'              spatial = TRUE,
##'              lrtest = FALSE,
##'              maxit = 20)
##'  }

asr = function(prep.out, fixed, random = ~1, spatial = TRUE, cor = TRUE, lrtest = FALSE,...) {
  
  stopifnot("'prep.out' must be an object of class 'comprepfor' or 'comprepcrop'" = class(prep.out) %in% c("comprepfor", "comprepcrop"))
  
  requireNamespace('asreml')
  
  control = attr(prep.out, 'control')
  
  prep.out <<- prep.out
  control <<- control
  
  input = prep.out$data
  
  if(inherits(prep.out, "comprepfor") && control[,6] > 0) ### Multi-areas ------------
  {
    input = input[order(input[,names(control)[6]],
                        input[,names(control)[3]], 
                        input[,names(control)[4]]),]
    
    if(spatial){
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
                               '~dsum(~ar1(', colnames(control)[3],'):ar1(',
                               colnames(control)[4],')|',
                               colnames(control)[6],')')
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
                               '~dsum(~ar1(', colnames(control)[3],'):ar1(',
                               colnames(control)[4],')|',
                               colnames(control)[6],')')
                             ),
                             data = input, ...)
      }
      
      
      if(lrtest & scm$converge){
        message("====> Starting likelihood ratio tests")
        
        bench = suppressWarnings({
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2],
                                        '+ grp(g1)'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~dsum(~ar1(', colnames(control)[3],'):ar1(',
                           colnames(control)[4],')|',
                           colnames(control)[6],')')
                         ),
                         data = input, ...)
        })
        
        
        lr.d = suppressWarnings({asreml::lrt(
          bench,
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+ grp(g1)'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~dsum(~ar1(', colnames(control)[3],'):ar1(',
                           colnames(control)[4],')|',
                           colnames(control)[6],')')
                         ),
                         data = input, ...)
        )
        })
        
        lr.i = suppressWarnings({asreml::lrt(
          bench,
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2]))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~dsum(~ar1(', colnames(control)[3],'):ar1(',
                           colnames(control)[4],')|',
                           colnames(control)[6],')')
                         ),
                         data = input, ...)
        )
        })
        
        lr.result = cbind(effect = c("DGE", 'IGE'), 
                          rbind(lr.d[,-1],lr.i[,-1]))
        rownames(lr.result) = NULL
        scm$lrt = lr.result
        message("====> LRT results:")
        print(lr.result)
      }
      
      
    }else{
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
                               '~dsum(~id(units)|',
                               colnames(control)[6],')')
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
                               '~dsum(~id(units)|',
                               colnames(control)[6],')')
                             ),
                             data = input, ...)
      }
      if(lrtest & scm$converge){
        message("====> Starting likelihood ratio tests")
        
        bench = suppressWarnings({
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2],
                                        '+ grp(g1)'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~dsum(~id(units)|',
                           colnames(control)[6],')')
                         ),
                         data = input, ...)
        })
        
        
        lr.d = suppressWarnings({asreml::lrt(
          bench,
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+ grp(g1)'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~dsum(~id(units)|',
                           colnames(control)[6],')')
                         ),
                         data = input, ...)
        )
        })
        
        lr.i = suppressWarnings({asreml::lrt(
          bench,
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2]))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~dsum(~id(units)|',
                           colnames(control)[6],')')
                         ),
                         data = input, ...)
        )
        })
        
        lr.result = cbind(effect = c("DGE", 'IGE'), 
                          rbind(lr.d[,-1],lr.i[,-1]))
        rownames(lr.result) = NULL
        scm$lrt = lr.result
        message("====> LRT results:")
        print(lr.result)
      }
    }

  } else ### Single area -----------
  {
    input = input[order(input[,names(control)[3]], 
                        input[,names(control)[4]]),]
    
    if(spatial){
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
                               '~ar1(', colnames(control)[3],'):ar1(',
                               colnames(control)[4],')')
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
                               '~ar1(', colnames(control)[3],'):ar1(',
                               colnames(control)[4],')')
                             ),
                             data = input, ...)
      }
      
      if(lrtest & scm$converge){
        message("====> Starting likelihood ratio tests")
        
        bench = suppressWarnings({
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2],
                                        '+ grp(g1)'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1(', colnames(control)[3],'):ar1(',
                           colnames(control)[4],')')
                         ),
                         data = input, ...)
        })
        
        
        lr.d = suppressWarnings({asreml::lrt(
          bench,
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+ grp(g1)'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1(', colnames(control)[3],'):ar1(',
                           colnames(control)[4],')')
                         ),
                         data = input, ...)
        )
        })
        
        lr.i = suppressWarnings({asreml::lrt(
          bench,
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2]))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1(', colnames(control)[3],'):ar1(',
                           colnames(control)[4],')')
                         ),
                         data = input, ...)
        )
        })
        
        lr.result = cbind(effect = c("DGE", 'IGE'), 
                          rbind(lr.d[,-1],lr.i[,-1]))
        rownames(lr.result) = NULL
        scm$lrt = lr.result
        message("====> LRT results:")
        print(lr.result)
      }
      

    }else{
      if(cor){
        scm = asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+ str(~', colnames(control)[2],
                                            '+ grp(g1), ~corh(2):id(',
                                            colnames(control)[2],'))'))
                             ),
                             group = list(g1 = 1:control[,2]),
                             data = input, ...)
      }else{
        scm = asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+', colnames(control)[2],
                                            '+ grp(g1)'))
                             ),
                             group = list(g1 = 1:control[,2]),
                             data = input, ...)
      }
      if(lrtest & scm$converge){
        message("====> Starting likelihood ratio tests")
        
        bench = suppressWarnings({
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2],
                                        '+ grp(g1)'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         data = input, ...)
        })
        
        
        lr.d = suppressWarnings({asreml::lrt(
          bench,
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+ grp(g1)'))
                         ),
                         group = list(g1 = 1:control[,2]),
                         data = input, ...)
        )
        })
        
        lr.i = suppressWarnings({asreml::lrt(
          bench,
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2]))
                         ),
                         group = list(g1 = 1:control[,2]),
                         data = input, ...)
        )
        })
        
        lr.result = cbind(effect = c("DGE", 'IGE'), 
                          rbind(lr.d[,-1],lr.i[,-1]))
        rownames(lr.result) = NULL
        scm$lrt = lr.result
        message("====> LRT results:")
        print(lr.result)
      }
    }
  }
  
  remove(prep.out, envir = .GlobalEnv)
  remove(control, envir = .GlobalEnv)
  class(scm) = c("compmod", class(scm))
  
  return(scm)
}





# Multi-environment spatial competition model ----------------------------

##' @title Fit a multi-age genetic competition model
##' 
##' @description
##' This function fits multi-age genetic competition models using [asreml::asreml()]. Currently
##' only models for forest breeding are implemented.
##' 
##' @param prep.out A `comprepfor` object.
##' @inheritParams asreml::asreml 
##' @param random The function is set to fit the genetic competition models in the 
##' `random` term. To consider further random effects, `random` receives a formula
##'  object specifying them. Otherwise, `random = ~1` (default). This argument has 
##'  the same general characteristics as fixed, but there can be no left side to 
##'  the ~ operator. Variance structures imposed on random terms are specified using 
##'  special model functions. See [asreml::asreml()] for details.
##' @param spatial A logical value. If `TRUE` (default), fits a spatial-genetic competition model 
##' @param cor A logical value. If `TRUE` (default), fits a model considering the correlation 
##' between direct and indirect genetic effects.
##' @param lrtest A logical value. If `TRUE`, performs a likelihood ratio test to verify
##' the significance of the direct and indirect genetic effects. Defaults to `FALSE`.
##' @inheritDotParams asreml::asreml -fixed -random -data -residual -group
##' 
##' @return An object of class `asreml` containing the results of the fitted linear model.
##' Instances of generic methods such as `plot()`, `predict()` and `summary()` return
##' various derived results of the fit. The method `resid()`, `coef()` and `fitted()` 
##' extract some of its components. See [asreml::asreml.object()] for details of the 
##' components of the returned list.
##'
##' @details
##' A general multi-age genetic competition linear mixed model can be represented by: 
##' 
##' \deqn{\mathbf{y} = \mathbf{X} \boldsymbol{\beta} + \mathbf{Z}_g \mathbf{g} + \mathbf{Z}_c \mathbf{c} + \mathbf{Z}_p \mathbf{p}  + \boldsymbol{\xi}}
##' 
##' where \eqn{\mathbf{y}} is the vector of phenotypic records, \eqn{\boldsymbol{\beta}}
##' is the vector of fixed effects, \eqn{\mathbf{g}} is the vector of direct genetic 
##' effects (DGE) within ages, \eqn{\mathbf{c}} is the vector of indirect genetic effects (IGE) within ages,
##' \eqn{\mathbf{p}} is the vector of other random effects, and \eqn{\boldsymbol{\varepsilon}} is 
##' the vector of errors. \eqn{\mathbf{X}} is the incidence matrix of the fixed effects, \eqn{\mathbf{Z}_g} 
##' is the DGE incidence matrix, \eqn{\mathbf{Z}_c} is the IGE incidence matrix , and \eqn{\mathbf{Z}_p} is the design matrix of 
##' other random effects. The dimensions of \eqn{\mathbf{Z}_c} are the same 
##' as \eqn{\mathbf{Z}_g}. \eqn{\mathbf{g}} and \eqn{\mathbf{c}} follow a compound symmetry structure, with 
##' a explicit variance for the main effects, and another for the genotype-by-environment interaction. 
##' The residual variance is assumed to be heterogeneous (one variance per age). 
##' If `spatial = TRUE`, the errors will be adjusted using a first-order autoregressive 
##' structure, with autocorrelation coefficients particularized per age. If `cor = TRUE`, the function will fit a model in which 
##' the main effects of \eqn{\mathbf{g}} and \eqn{\mathbf{c}} are correlated outcomes of 
##' the main genotypic effects decomposition. They both follow a Gaussian distribution, 
##' with mean centred in zero, and covariance given by:
##'  
##'  \deqn{\mathbf{\Sigma_g} =  \begin{bmatrix}\sigma_{\text{g}}^2 & \sigma_{\text{gc}}\\\sigma_{\text{gc}} & \sigma_{\text{c}}^2\\\end{bmatrix}\otimes {{\mathbf I_V}}}
##'  
##'  where \eqn{\sigma_{\text{g}}^2} is the DGE variance, \eqn{\sigma_{\text{c}}^2} is the IGE variance, 
##'  and \eqn{\sigma_{\text{gc}}} is the covariance between DGE and IGE.
##'  
##'  The likelihood ratio test is performed using a model without the correlation between DGE and IGE.
##'  
##'
##' @seealso  [gencomp::prepfor], [asreml::asreml.options], [asreml::asreml.object], [asreml::family_dist]
##'
##' @import asreml
##' @importFrom stats as.formula
##' 
##' @export
##' 
##' @examples
##' \donttest{
##'  library(gencomp)
##'  comp_mat = prepfor(data = euca, gen = 'clone', area = 'area',
##'                    ind = 'tree', age = 'age', row = 'row', col = 'col',
##'                    dist.col = 3, dist.row = 2, trait = 'MAI', method = 'SK',
##'                    n.dec = 3, verbose = FALSE, effs = c("block"))
##'  model = asr_ma(prep.out = comp_mat,
##'                 fixed = MAI~ age, 
##'                 random = ~ block:age, 
##'                 lrtest = TRUE, 
##'                 spatial = TRUE, 
##'                 cor = TRUE, 
##'                 maxit = 20)
##'  }

asr_ma = function(prep.out, fixed, random = ~1, spatial = TRUE, cor = TRUE,
                  lrtest = FALSE,...) {
  
  stopifnot("'prep.out' must be an object of class 'comprepfor' or 'comprepcrop'" = class(prep.out) %in% c("comprepfor", "comprepcrop"))
  stopifnot("Currently, only models for 'comprepfor' objects are implemented" = class(prep.out) == 'comprepfor')
  requireNamespace('asreml')
  
  control = attr(prep.out, 'control')
  
  prep.out <<- prep.out
  control <<- control
  
  input = prep.out$data
  
  if(control[,6] > 0) ### Multi-areas ------------
  {
    input[,"aux"] = factor(paste(paste(names(control)[6], input[,names(control)[6]], sep = '_'),
                                 paste(names(control)[5], input[,names(control)[5]], sep = '_'),sep = ':'))
    input = input[order(input[,'aux'], input[,names(control)[3]], input[,names(control)[4]]),]
    
    if(spatial){
        
        if(cor){
          scm = asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+ str(~', colnames(control)[2],
                                              '+ grp(g1), ~corh(2):id(',
                                              colnames(control)[2],')) + str(~',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5], ', ~diag(2):id(',
                                              colnames(control)[5], '):id(',
                                              colnames(control)[2], '))'))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(paste0(
                                 '~dsum(~ar1(',
                                 colnames(control)[3],'):ar1(',
                                 colnames(control)[4],')|aux)')
                               ),
                               data = input, ...)
        }else{
          scm = asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+', colnames(control)[2],
                                              '+ grp(g1) +',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(paste0(
                                 '~dsum(~ar1(',
                                 colnames(control)[3],'):ar1(',
                                 colnames(control)[4],')|aux)')
                               ),
                               data = input, ...)
          
        }
        
      if(lrtest & scm$converge){
        message("====> Starting likelihood ratio tests")
        
        bench = suppressWarnings({
          asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(control)[2],
                                        '+ grp(g1) +',
                                        colnames(control)[2],':',
                                        colnames(control)[5],'+grp(g1):',
                                        colnames(control)[5]))
                         ),
                         group = list(g1 = 1:control[,2]),
                         residual = stats::as.formula(paste0(
                           '~dsum(~ar1(',
                           colnames(control)[3],'):ar1(',
                           colnames(control)[4],')|aux)')
                         ),
                         data = input, ...)
        }) 
        
        lr.dint = suppressWarnings({
          lrt(bench, 
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+', colnames(control)[2],
                                            '+ grp(g1) +','grp(g1):',
                                            colnames(control)[5]))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~dsum(~ar1(',
                               colnames(control)[3],'):ar1(',
                               colnames(control)[4],')|aux)')
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
                                            colnames(control)[5]))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~dsum(~ar1(',
                               colnames(control)[3],'):ar1(',
                               colnames(control)[4],')|aux)')
                             ),
                             data = input, ...))
        }) 
        
        lr.d = suppressWarnings({
          lrt(bench, 
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+ grp(g1) +',
                                            colnames(control)[2],':',
                                            colnames(control)[5],'+grp(g1):',
                                            colnames(control)[5]))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~dsum(~ar1(',
                               colnames(control)[3],'):ar1(',
                               colnames(control)[4],')|aux)')
                             ),
                             data = input, ...))
        })
        
        lr.i = suppressWarnings({
          lrt(bench, 
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+', colnames(control)[2],
                                            '+',
                                            colnames(control)[2],':',
                                            colnames(control)[5],'+grp(g1):',
                                            colnames(control)[5]))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~dsum(~ar1(',
                               colnames(control)[3],'):ar1(',
                               colnames(control)[4],')|aux)')
                             ),
                             data = input, ...))
        }) 
        
        lr.result = cbind(effect = c("DGE", 'IGE', 'DGE:age', 'IGE:age'), 
                          rbind(lr.d[,-1], lr.i[,-1], lr.dint[,-1], lr.iint[,-1]))
        rownames(lr.result) = NULL
        
        message("====> LRT results:")
        print(lr.result)
        scm$lrt = lr.result
      }
      } else {
        if(cor){
          scm = asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+ str(~', colnames(control)[2],
                                              '+ grp(g1), ~corh(2):id(',
                                              colnames(control)[2],')) + str(~',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5], ', ~diag(2):id(',
                                              colnames(control)[5], '):id(',
                                              colnames(control)[2], '))'))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = ~dsum(~id(units)|aux),
                               data = input, ...)
        }else{
          scm = asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+', colnames(control)[2],
                                              '+ grp(g1) +',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = ~dsum(~id(units)|aux),
                               data = input, ...)
          
        }
        if(lrtest & scm$converge){
          message("====> Starting likelihood ratio tests")
          
          bench = suppressWarnings({
            asreml::asreml(fixed = fixed,
                           random = stats::as.formula(
                             paste(paste(as.character(random), collapse = ""), 
                                   paste0('+', colnames(control)[2],
                                          '+ grp(g1) +',
                                          colnames(control)[2],':',
                                          colnames(control)[5],'+grp(g1):',
                                          colnames(control)[5]))
                           ),
                           group = list(g1 = 1:control[,2]),
                           residual = ~dsum(~id(units)|aux),
                           data = input, ...)
          }) 
          
          lr.dint = suppressWarnings({
            lrt(bench, 
                asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+', colnames(control)[2],
                                              '+ grp(g1) +','grp(g1):',
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = ~dsum(~id(units)|aux),
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
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = ~dsum(~id(units)|aux),
                               data = input, ...))
          }) 
          
          lr.d = suppressWarnings({
            lrt(bench, 
                asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+ grp(g1) +',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = ~dsum(~id(units)|aux),
                               data = input, ...))
          })
          
          lr.i = suppressWarnings({
            lrt(bench, 
                asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+', colnames(control)[2],
                                              '+',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = ~dsum(~id(units)|aux),
                               data = input, ...))
          }) 
          
          lr.result = cbind(effect = c("DGE", 'IGE', 'DGE:age', 'IGE:age'), 
                            rbind(lr.d[,-1], lr.i[,-1], lr.dint[,-1], lr.iint[,-1]))
          rownames(lr.result) = NULL
          
          message("====> LRT results:")
          print(lr.result)
          scm$lrt = lr.result
        }
      }
      
  } else ### Single area -----------
  {
      input = input[order(input[,names(control)[5]],
                          input[,names(control)[3]], 
                          input[,names(control)[4]]),]
      
      if(spatial){
        if(cor){
          scm = asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+ str(~', colnames(control)[2],
                                              '+ grp(g1), ~corh(2):id(',
                                              colnames(control)[2],')) + str(~',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5], ', ~diag(2):id(',
                                              colnames(control)[5], '):id(',
                                              colnames(control)[2], '))'))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(
                                 paste0(
                                   "~dsum(~ar1(",colnames(control)[3],
                                   "):ar1(",colnames(control)[4],")|",
                                   colnames(control)[5],")"
                                 )
                               ),
                               data = input, ...)
          
        }else{
          scm = asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+', colnames(control)[2],
                                              '+ grp(g1) +',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(
                                 paste0(
                                   "~dsum(~ar1(",colnames(control)[3],
                                   "):ar1(",colnames(control)[4],")|",
                                   colnames(control)[5],")"
                                 )
                               ),
                               data = input, ...)
        }
        
        if(lrtest & scm$converge){
          message("====> Starting likelihood ratio tests")
          
          bench = suppressWarnings({
            asreml::asreml(fixed = fixed,
                           random = stats::as.formula(
                             paste(paste(as.character(random), collapse = ""), 
                                   paste0('+', colnames(control)[2],
                                          '+ grp(g1) +',
                                          colnames(control)[2],':',
                                          colnames(control)[5],'+grp(g1):',
                                          colnames(control)[5]))
                           ),
                           group = list(g1 = 1:control[,2]),
                           residual = stats::as.formula(
                             paste0(
                               "~dsum(~ar1(",colnames(control)[3],
                               "):ar1(",colnames(control)[4],")|",
                               colnames(control)[5],")"
                             )
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
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(
                                 paste0(
                                   "~dsum(~ar1(",colnames(control)[3],
                                   "):ar1(",colnames(control)[4],")|",
                                   colnames(control)[5],")"
                                 )
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
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(
                                 paste0(
                                   "~dsum(~ar1(",colnames(control)[3],
                                   "):ar1(",colnames(control)[4],")|",
                                   colnames(control)[5],")"
                                 )
                               ),
                               data = input, ...))
          }) 
          
          lr.d = suppressWarnings({
            lrt(bench, 
                asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+ grp(g1) +',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(
                                 paste0(
                                   "~dsum(~ar1(",colnames(control)[3],
                                   "):ar1(",colnames(control)[4],")|",
                                   colnames(control)[5],")"
                                 )
                               ),
                               data = input, ...))
          })
          
          lr.i = suppressWarnings({
            lrt(bench, 
                asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+', colnames(control)[2],
                                              '+',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(
                                 paste0(
                                   "~dsum(~ar1(",colnames(control)[3],
                                   "):ar1(",colnames(control)[4],")|",
                                   colnames(control)[5],")"
                                 )
                               ),
                               data = input, ...))
          }) 
          
          lr.result = cbind(effect = c("DGE", 'IGE', 'DGE:age', 'IGE:age'), 
                            rbind(lr.d[,-1], lr.i[,-1], lr.dint[,-1], lr.iint[,-1]))
          rownames(lr.result) = NULL
          
          message("====> LRT results:")
          print(lr.result)
          scm$lrt = lr.result
        }
        
        
      } else {
        if(cor){
          scm = asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+ str(~', colnames(control)[2],
                                              '+ grp(g1), ~corh(2):id(',
                                              colnames(control)[2],')) + str(~',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5], ', ~diag(2):id(',
                                              colnames(control)[5], '):id(',
                                              colnames(control)[2], '))'))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(
                                 paste0(
                                   "~dsum(~id(units)|",
                                   colnames(control)[5],")"
                                 )
                               ),
                               data = input, ...)
          
        }else{
          scm = asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+', colnames(control)[2],
                                              '+ grp(g1) +',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(
                                 paste0(
                                   "~dsum(~id(units)|",
                                   colnames(control)[5],")"
                                 )
                               ),
                               data = input, ...)
        }
        
        if(lrtest & scm$converge){
          message("====> Starting likelihood ratio tests")
          
          bench = suppressWarnings({
            asreml::asreml(fixed = fixed,
                           random = stats::as.formula(
                             paste(paste(as.character(random), collapse = ""), 
                                   paste0('+', colnames(control)[2],
                                          '+ grp(g1) +',
                                          colnames(control)[2],':',
                                          colnames(control)[5],'+grp(g1):',
                                          colnames(control)[5]))
                           ),
                           group = list(g1 = 1:control[,2]),
                           residual = stats::as.formula(
                             paste0(
                               "~dsum(~id(units)|",
                               colnames(control)[5],")"
                             )
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
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(
                                 paste0(
                                   "~dsum(~id(units)|",
                                   colnames(control)[5],")"
                                 )
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
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(
                                 paste0(
                                   "~dsum(~id(units)|",
                                   colnames(control)[5],")"
                                 )
                               ),
                               data = input, ...))
          }) 
          
          lr.d = suppressWarnings({
            lrt(bench, 
                asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+ grp(g1) +',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(
                                 paste0(
                                   "~dsum(~id(units)|",
                                   colnames(control)[5],")"
                                 )
                               ),
                               data = input, ...))
          })
          
          lr.i = suppressWarnings({
            lrt(bench, 
                asreml::asreml(fixed = fixed,
                               random = stats::as.formula(
                                 paste(paste(as.character(random), collapse = ""), 
                                       paste0('+', colnames(control)[2],
                                              '+',
                                              colnames(control)[2],':',
                                              colnames(control)[5],'+grp(g1):',
                                              colnames(control)[5]))
                               ),
                               group = list(g1 = 1:control[,2]),
                               residual = stats::as.formula(
                                 paste0(
                                   "~dsum(~id(units)|",
                                   colnames(control)[5],")"
                                 )
                               ),
                               data = input, ...))
          }) 
          
          lr.result = cbind(effect = c("DGE", 'IGE', 'DGE:age', 'IGE:age'), 
                            rbind(lr.d[,-1], lr.i[,-1], lr.dint[,-1], lr.iint[,-1]))
          rownames(lr.result) = NULL
          
          message("====> LRT results:")
          print(lr.result)
          scm$lrt = lr.result
        }
        
      }
  }
  
  remove(prep.out, envir = .GlobalEnv)
  remove(control, envir = .GlobalEnv)
  
  class(scm) = c("compmod", class(scm))
  
  return(scm)

}