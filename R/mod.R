##' @title Fit a genetic-spatial competition model
##' 
##' @description
##' This function fits genetic-spatial competition models using ASReml-R
##' 
##' @param prep.out A `comp.prep` object.
##' @inheritParams asreml::asreml 
##' @param random The function is set to fit the genetic competition models in the 
##' `random` term. To consider further random effects, `random` receives a formula
##'  object specifying them. Otherwise, `random = ~1` (default). This argument has 
##'  the same general characteristics as fixed, but there can be no left side to 
##'  the ~ operator. Variance structures imposed on random terms are specified using 
##'  special model functions.
##' @param cor A logical value. If `TRUE` (default), fits a model considering the correlation 
##' between direct and indirect genetic effects.
##' @param maxit An integer. The maximum number of iterations. Defaults to 50.
##' @inheritDotParams asreml::asreml -fixed -random -data
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
##' using [competition::comp.prep]), and \eqn{\mathbf{Z}_p} is the design matrix of 
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
##'  If the dataset has multiple age, the model described aboved is expanded with two more effects:
##'  the DGE \eqn{\times} age interaction, and the IGE \eqn{\times} age interaction.
##'  
##'
##' @seealso  [competition::comp.prep], [asreml::asreml.options()], [asreml::asreml.object()], [asreml::family_dist()]
##'
##' @import asreml
##' @importFrom stats as.formula
##' 
##' @export
##' 
##' @examples
##' \donttest{
##'  comp_mat = comp.prep(data = data, gen = 'clone', repl = 'block', area = 'area', 
##'                       ind = 'tree', age = 'age', row = 'row', col = 'col', 
##'                       dist.col = 3, dist.row = 2.5, trait = 'mai', method = 'SK',
##'                       n.dec = 3, verbose = TRUE)
##'  
##'  model = comp.asr(prep.out = comp_mat, 
##'                   fixed = mai~ age, 
##'                   random = ~ block:age, 
##'                   cor = T, 
##'                   maxit = 50)
##'  }

comp.asr = function(prep.out, fixed, random = ~1, cor = TRUE, maxit = 50, ...) {
  
  requireNamespace('asreml')
  
  prep.out <<- prep.out
  
  input = prep.out$data
  
  if(prep.out$control[,7] > 0) ### Multi-areas ------------
    {
      if(prep.out$control[,6] > 0) ### Multi-ages --------
      {
        if(cor){
          input = input[order(input$area, input$age, input$row, input$col),]
          
          scm = asreml::asreml(fixed = fixed,
                       random = stats::as.formula(
                         paste(paste(as.character(random), collapse = ""), 
                               paste0('+ str(~', colnames(prep.out$control)[2],
                                      '+ grp(g1), ~corh(2):id(',
                                      colnames(prep.out$control)[2],')) + str(~',
                                      colnames(prep.out$control)[2],':',
                                      colnames(prep.out$control)[6],'+grp(g1):',
                                      colnames(prep.out$control)[6], ', ~diag(2):id(',
                                      colnames(prep.out$control)[6], '):id(',
                                      colnames(prep.out$control)[2], '))'))
                       ),
                       group = list(g1 = 1:prep.out$control[,2]),
                       residual = stats::as.formula(paste0(
                         '~dsum(~ar1h(', colnames(prep.out$control)[6],'):ar1(',
                         colnames(prep.out$control)[4],'):ar1(',
                         colnames(prep.out$control)[5],')|',
                         colnames(prep.out$control)[7],')')
                       ),
                       maxit = maxit,
                       data = input, ...)

        }else{
          input = input[order(input$area, input$age, input$row, input$col),]
          
          scm = asreml::asreml(fixed = fixed,
                       random = stats::as.formula(
                         paste(paste(as.character(random), collapse = ""), 
                               paste0('+', colnames(prep.out$control)[2],
                                      '+ grp(g1) +',
                                      colnames(prep.out$control)[2],':',
                                      colnames(prep.out$control)[6],'+grp(g1):',
                                      colnames(prep.out$control)[6]))
                       ),
                       group = list(g1 = 1:prep.out$control[,2]),
                       residual = stats::as.formula(paste0(
                         '~dsum(~ar1h(', colnames(prep.out$control)[6],'):ar1(',
                         colnames(prep.out$control)[4],'):ar1(',
                         colnames(prep.out$control)[5],')|',
                         colnames(prep.out$control)[7],')')
                       ),
                       maxit = maxit,
                       data = input, ...)

        }
      }else ### Single age ---------------
      {
        if(cor){
          scm = asreml::asreml(fixed = fixed,
                       random = stats::as.formula(
                         paste(paste(as.character(random), collapse = ""), 
                               paste0('+ str(~', colnames(prep.out$control)[2],
                                      '+ grp(g1), ~corh(2):id(',
                                      colnames(prep.out$control)[2],'))'))
                       ),
                       group = list(g1 = 1:prep.out$control[,2]),
                       residual = stats::as.formula(paste0(
                         '~dsum(~ar1(', colnames(prep.out$control)[4],'):ar1(',
                         colnames(prep.out$control)[5],')|',
                         colnames(prep.out$control)[7],')')
                       ),
                       maxit = maxit,
                       data = input, ...)
        }else{
          scm = asreml::asreml(fixed = fixed,
                       random = stats::as.formula(
                         paste(paste(as.character(random), collapse = ""), 
                               paste0('+', colnames(prep.out$control)[2],
                                      '+ grp(g1)'))
                       ),
                       group = list(g1 = 1:prep.out$control[,2]),
                       residual = stats::as.formula(paste0(
                         '~dsum(~ar1(', colnames(prep.out$control)[4],'):ar1(',
                         colnames(prep.out$control)[5],')|',
                         colnames(prep.out$control)[7],')')
                       ),
                       maxit = maxit,
                       data = input, ...)
        }
      } 
    }else ### Single area -----------
      {
        if(prep.out$control[,6] > 0) ### Multi-ages --------
        {
          if(cor){
            input = input[order(input$age, input$row, input$col),]
            scm = asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+ str(~', colnames(prep.out$control)[2],
                                        '+ grp(g1), ~corh(2):id(',
                                        colnames(prep.out$control)[2],')) + str(~',
                                        colnames(prep.out$control)[2],':',
                                        colnames(prep.out$control)[6],'+grp(g1):',
                                        colnames(prep.out$control)[6], ', ~corh(2):id(',
                                        colnames(prep.out$control)[6], '):id(',
                                        colnames(prep.out$control)[2], '))'))
                         ),
                         group = list(g1 = 1:prep.out$control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1h(', colnames(prep.out$control)[6],'):ar1(',
                           colnames(prep.out$control)[4],'):ar1(',
                           colnames(prep.out$control)[5],')')
                         ),
                         maxit = maxit,
                         data = input, ...)
            
          }else{
            input = input[order(input$age, input$row, input$col),]
            scm = asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(prep.out$control)[2],
                                        '+ grp(g1) +',
                                        colnames(prep.out$control)[2],':',
                                        colnames(prep.out$control)[6],'+grp(g1):',
                                        colnames(prep.out$control)[6]))
                         ),
                         group = list(g1 = 1:prep.out$control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1h(', colnames(prep.out$control)[6],'):ar1(',
                           colnames(prep.out$control)[4],'):ar1(',
                           colnames(prep.out$control)[5],')')
                         ),
                         maxit = maxit,
                         data = input, ...)
          }
        }else ### Single age ---------------
        {
          if(cor){
            scm = asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+ str(~', colnames(prep.out$control)[2],
                                        '+ grp(g1), ~corh(2):id(',
                                        colnames(prep.out$control)[2],'))'))
                         ),
                         group = list(g1 = 1:prep.out$control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1(', colnames(prep.out$control)[4],'):ar1(',
                           colnames(prep.out$control)[5],')')
                         ),
                         maxit = maxit,
                         data = input, ...)
          }else{
            scm = asreml::asreml(fixed = fixed,
                         random = stats::as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(prep.out$control)[2],
                                        '+ grp(g1)'))
                         ),
                         group = list(g1 = 1:prep.out$control[,2]),
                         residual = stats::as.formula(paste0(
                           '~ar1(', colnames(prep.out$control)[4],'):ar1(',
                           colnames(prep.out$control)[5],')')
                         ),
                         maxit = maxit,
                         data = input, ...)
          }
        }
  }
  
  if(!scm$converge) warning("The model did not converge")
  
  remove(prep.out, envir = .GlobalEnv)
  
  return(scm)
}


