# Update asreml model ---------------

##' @title Update the asreml object
##' 
##' @description
##' Updates the `asreml` object until full convergence
##' 
##' @param model An object of class `asreml`.
##' @return An object of class `asreml` with all parameters converged.
##' @seealso [gencomp::asr()], [asreml::asreml.options], [asreml::asreml.object]
##' 
##' @importFrom asreml update.asreml
##' @importFrom stats na.exclude
##' @keywords internal

.up.mod = function(model) {
  if (model$converge) {
    repeat {
      if (any(stats::na.exclude(model$vparameters.pc) >= 1)) {
        model = suppressWarnings(asreml::update.asreml(model))
      } else{
        message("The model was updated and reached full convergence")
        break
      }
    }
  } else message("The model failed to converge")
  return(model)
}

# Genetic competition model ----------------------------

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
##' @param K A sparse inverse variance matrix held in three-column coordinate form in row-major order. 
##' For backwards compatibility, a three-column data frame is also accepted. 
##' In either case, the source must have a `rowNames` attribute. 
##' This object must have class `ginv` from a call to [asreml::ainverse()], 
##' or have attribute `INVERSE` set to `TRUE`. If the attribute `INVERSE` is not set then FALSE is assumed. 
##' We recommend using the `ASRgenomics` package to set the genomic relationship matrix (see [ASRgenomics::G.matrix()] and 
##' [ASRgenomics::G.inverse()]). If there is no relationship information, `K = NULL` (default)
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
##'  \eqn{\mathbf{g}} and \eqn{\mathbf{c}} are correlated outcomes of the genetic effects'
##'  decomposition. They both follow a Gaussian distribution, with mean centred in zero, and
##'  covariance given by:
##'  
##'  \deqn{\mathbf{\Sigma_g} = \begin{bmatrix}\sigma_{\text{g}}^2 & \sigma_{\text{gc}}\\\sigma_{\text{gc}} & \sigma_{\text{c}}^2\\\end{bmatrix}\otimes {{\mathbf K}}}
##'  
##'  where \eqn{\sigma_{\text{g}}^2} is the DGE variance, \eqn{\sigma_{\text{c}}^2} is the IGE variance, 
##'  \eqn{\sigma_{\text{gc}}} is the covariance between DGE and IGE, and \eqn{K} is a \eqn{V \times V} relationship matrix.
##'  
##'  If `prep.out` is a `comprepfor` object, and the argument `age` is different than `NULL` in 
##'  the [gencomp::prepfor] function, a multi-age genetic competition model is fitted. In this case, 
##'  the vectors \eqn{\mathbf{g}} and \eqn{\mathbf{c}} follow a compound symmetry structure, with 
##'  a explicit variance for the main effects, and another for the genotype-by-environment interaction. 
##'  The residual variance is assumed to be heterogeneous (one variance per age). 
##'  If `spatial = TRUE`, the errors will be adjusted using a first-order autoregressive 
##'  structure, with autocorrelation coefficients particularized per age. If `cor = TRUE`, the function will fit a model in which 
##'  the main effects of \eqn{\mathbf{g}} and \eqn{\mathbf{c}} are correlated outcomes of 
##'  the main genetic effects decomposition. 
##'  
##'  The likelihood ratio test is performed comparing models excluding the covariance component.
##'  
##'
##' @seealso  [gencomp::prepfor], [gencomp::prepcrop], [asreml::asreml.options], [asreml::asreml.object], [asreml::ainverse]
##'
##' @import asreml
##' @importFrom stats as.formula
##' 
##' @export
##' 
##' @examples
##' \donttest{
##'  library(gencomp)
##'  
##'  # Genetic competition model using a crop breeding example
#'   comps = prepcrop(data = potato, gen = "gen", row = "row", col = "col",
#'                    plt = NULL, effs = c("rep", 'matur'), trait = "yield",
#'                    direction = "row", verbose = TRUE)
#'  
#'   mod = asr(prep.out = comps, 
#'             fixed = yield ~ matur + rep,
#'             random = ~1, 
#'             spatial = FALSE, 
#'             cor = TRUE, 
#'             lrtest = TRUE, 
#'             K = NULL)
#'   
#'   # Genetic competition model using a forest tree breeding example
#'   comp_mat = prepfor(data = euca, gen = 'clone', area = 'area',
##'                    plt = 'tree', age = 'age', row = 'row', col = 'col',
##'                    dist.col = 3, dist.row = 2, trait = 'MAI', method = 'SK',
##'                    n.dec = 3, verbose = FALSE, effs = c("block"))
##'  model = asr(prep.out = comp_mat,
##'              fixed = MAI ~ age, 
##'              random = ~ block:age, 
##'              lrtest = TRUE, 
##'              spatial = TRUE, 
##'              cor = TRUE, 
##'              K = NULL,
##'              maxit = 20)
##'  }

asr = function(prep.out, fixed, random = ~1, spatial = TRUE, cor = TRUE, 
               lrtest = FALSE, K = NULL, ...) {
  
  stopifnot("'prep.out' must be an object of class 'comprepfor' or 'comprepcrop'" = class(prep.out) %in% c("comprepfor", "comprepcrop"))

  requireNamespace('asreml')
  
  input = prep.out$data
  control = attr(prep.out, 'control')
  fixed <<- fixed
  random <<- random
  
  ## Adjusting the data frame if there is relationship information
  if(!is.null(K)){
    stopifnot("'K' must have a rowNames attribute" = "rowNames" %in% names(attributes(K)))
    trats = attr(K, "rowNames")
    founders = sort(trats[which(!trats %in% input[,names(control)[2]])])
    founders_matrix = matrix(0, nrow = nrow(input), ncol = length(founders),
                             dimnames = list(NULL, founders))
    input = data.frame(founders_matrix, input)
    input = input[, c(trats, colnames(input)[which(!colnames(input) %in% trats)])]
    control[,2] = control[,2] + length(founders)
    K <<- K
  }
  control <<- control
  
  if(inherits(prep.out, "comprepfor") && control[5] > 0) multi.age = TRUE else multi.age = FALSE
  if(inherits(prep.out, "comprepfor") && multi.age) stopifnot("Currently, multi-age models are implemented only for 'comprepfor' objects" = class(prep.out) == 'comprepfor')
  
 
  if(multi.age){ 
    if(control[,6] > 1) ### Multi-areas ------------
    {
      input[,"aux"] = factor(paste(paste(names(control)[6], input[,names(control)[6]], sep = '_'),
                                   paste(names(control)[5], input[,names(control)[5]], sep = '_'),sep = ':'))
      input <<- input[order(input[,'aux'], input[,names(control)[3]], input[,names(control)[4]]),]
      
      if (spatial) {
        if (cor) {
          if (!is.null(K)) {
            scm = asreml::asreml(
              fixed = fixed,
              random = stats::as.formula(paste(
                paste(as.character(random), collapse = ""),
                paste0(
                  '+ str(~vm(',
                  colnames(control)[2],
                  ',K)+ grp(g1), ~corh(2):vm(',
                  colnames(control)[2],
                  ',K)) + str(~vm(',
                  colnames(control)[2],
                  ',K):',
                  colnames(control)[5],
                  '+grp(g1):',
                  colnames(control)[5],
                  ', ~diag(2):id(',
                  colnames(control)[5],
                  '):vm(',
                  colnames(control)[2],
                  ',K))'
                )
              )),
              group = list(g1 = 1:control[, 2]),
              residual = stats::as.formula(
                paste0(
                  '~dsum(~ar1(',
                  colnames(control)[3],
                  '):ar1(',
                  colnames(control)[4],
                  ')|aux)'
                )
              ),
              data = input,
              ...
            )
          } else{
            scm = asreml::asreml(
              fixed = fixed,
              random = stats::as.formula(paste(
                paste(as.character(random), collapse = ""),
                paste0(
                  '+ str(~',
                  colnames(control)[2],
                  '+ grp(g1), ~corh(2):id(',
                  colnames(control)[2],
                  ')) + str(~',
                  colnames(control)[2],
                  ':',
                  colnames(control)[5],
                  '+grp(g1):',
                  colnames(control)[5],
                  ', ~diag(2):id(',
                  colnames(control)[5],
                  '):id(',
                  colnames(control)[2],
                  '))'
                )
              )),
              group = list(g1 = 1:control[, 2]),
              residual = stats::as.formula(
                paste0(
                  '~dsum(~ar1(',
                  colnames(control)[3],
                  '):ar1(',
                  colnames(control)[4],
                  ')|aux)'
                )
              ),
              data = input,
              ...
            )
          }
        } else{
          if (!is.null(K)) {
            scm = asreml::asreml(
              fixed = fixed,
              random = stats::as.formula(paste(
                paste(as.character(random), collapse = ""),
                paste0(
                  '+ vm(',
                  colnames(control)[2],
                  ',K)+ grp(g1) + vm(',
                  colnames(control)[2],
                  ',K):',
                  colnames(control)[5],
                  '+grp(g1):',
                  colnames(control)[5]
                )
              )),
              group = list(g1 = 1:control[, 2]),
              residual = stats::as.formula(
                paste0(
                  '~dsum(~ar1(',
                  colnames(control)[3],
                  '):ar1(',
                  colnames(control)[4],
                  ')|aux)'
                )
              ),
              data = input,
              ...
            )
          } else{
            scm = asreml::asreml(
              fixed = fixed,
              random = stats::as.formula(paste(
                paste(as.character(random), collapse = ""),
                paste0(
                  '+',
                  colnames(control)[2],
                  '+ grp(g1) +',
                  colnames(control)[2],
                  ':',
                  colnames(control)[5],
                  '+grp(g1):',
                  colnames(control)[5]
                )
              )),
              group = list(g1 = 1:control[, 2]),
              residual = stats::as.formula(
                paste0(
                  '~dsum(~ar1(',
                  colnames(control)[3],
                  '):ar1(',
                  colnames(control)[4],
                  ')|aux)'
                )
              ),
              data = input,
              ...
            )
          }
        }
        
        if (any(na.exclude(scm$vparameters.pc > 1))) scm = .up.mod(scm)
        
        if(lrtest & scm$converge){
          message("====> Starting likelihood ratio tests")
          
          if (!is.null(K)){
            bench = suppressWarnings({
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+vm(', colnames(control)[2],
                                            ',K)+ grp(g1) +vm(',
                                            colnames(control)[2],',K):',
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
                                         paste0('+vm(', colnames(control)[2],
                                                ',K)+ grp(g1) +','grp(g1):',
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
                                         paste0('+vm(', colnames(control)[2],
                                                ',K)+ grp(g1) +',
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
                                         paste0('+ grp(g1)','+grp(g1):',
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
                                         paste0('+vm(', colnames(control)[2],
                                                ',K)+vm(',
                                                colnames(control)[2],',K):',
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
          }else{
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
                                         paste0('+ grp(g1)','+grp(g1):',
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
          }
          
          lr.result = cbind(effect = c("DGE", 'IGE', 'DGE:age', 'IGE:age'), 
                            rbind(lr.d[,-1], lr.i[,-1], lr.dint[,-1], lr.iint[,-1]))
          rownames(lr.result) = NULL
          
          message("====> LRT results:")
          print(lr.result)
          scm$lrt = lr.result
        }
      } else {
        if(cor){
          if (!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ str(~vm(', colnames(control)[2],
                                                ',K)+ grp(g1), ~corh(2):vm(',
                                                colnames(control)[2],',K)) + str(~vm(',
                                                colnames(control)[2],',K):',
                                                colnames(control)[5],'+grp(g1):',
                                                colnames(control)[5], ', ~diag(2):id(',
                                                colnames(control)[5], '):vm(',
                                                colnames(control)[2], ',K))'))
                                 ),
                                 group = list(g1 = 1:control[,2]),
                                 residual = ~dsum(~id(units)|aux),
                                 data = input, ...)
          }else{
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
          }
          
        }else{
          if (!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K)+ grp(g1) + vm(',
                                                colnames(control)[2],',K):',
                                                colnames(control)[5],' + grp(g1):',
                                                colnames(control)[5]))
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
        }
        
        if (any(na.exclude(scm$vparameters.pc > 1))) scm = .up.mod(scm)
        
        if(lrtest & scm$converge){
          message("====> Starting likelihood ratio tests")
          
          if (!is.null(K)){
            bench = suppressWarnings({
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+ vm(', colnames(control)[2],
                                            ',K)+ grp(g1) + vm(',
                                            colnames(control)[2],',K):',
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
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K) + grp(g1) +','grp(g1):',
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
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K) + grp(g1) +',
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
                                         paste0('+ grp(g1)','+ grp(g1):',
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
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K) + vm(',
                                                colnames(control)[2],',K):',
                                                colnames(control)[5]))
                                 ),
                                 group = list(g1 = 1:control[,2]),
                                 residual = ~dsum(~id(units)|aux),
                                 data = input, ...))
            }) 
          }else{
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
                                         paste0('+ grp(g1)','+ grp(g1):',
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
                                                colnames(control)[5]))
                                 ),
                                 group = list(g1 = 1:control[,2]),
                                 residual = ~dsum(~id(units)|aux),
                                 data = input, ...))
            }) 
          }
          
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
      input <<- input[order(input[,names(control)[5]],
                            input[,names(control)[3]], 
                            input[,names(control)[4]]),]
      
      if(spatial){
        if(cor){
          if (!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ str(~ vm(', colnames(control)[2],
                                                ',K)+ grp(g1), ~corh(2):vm(',
                                                colnames(control)[2],',K)) + str(~vm(',
                                                colnames(control)[2],',K):',
                                                colnames(control)[5],'+grp(g1):',
                                                colnames(control)[5], ', ~diag(2):id(',
                                                colnames(control)[5], '):vm(',
                                                colnames(control)[2], ',K))'))
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
          }
          
        }else{
          if (!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K) + grp(g1) + vm(',
                                                colnames(control)[2],',K):',
                                                colnames(control)[5],' + grp(g1):',
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
          
        }
        
        if (any(na.exclude(scm$vparameters.pc > 1))) scm = .up.mod(scm)
        
        if(lrtest & scm$converge){
          message("====> Starting likelihood ratio tests")
          
          if (!is.null(K)){
            bench = suppressWarnings({
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+ vm(', colnames(control)[2],
                                            ',K) + grp(g1) + vm(',
                                            colnames(control)[2],',K):',
                                            colnames(control)[5],' + grp(g1):',
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
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K) + grp(g1)','+grp(g1):',
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
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K) + grp(g1) + vm(',
                                                colnames(control)[2],',K):',
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
                                         paste0('+ grp(g1) ','+ grp(g1):',
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
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K) + vm(',
                                                colnames(control)[2],',K):',
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
          }else{
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
                                         paste0('+ grp(g1) +',' grp(g1):',
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
          }
          
          lr.result = cbind(effect = c("DGE", 'IGE', 'DGE:age', 'IGE:age'), 
                            rbind(lr.d[,-1], lr.i[,-1], lr.dint[,-1], lr.iint[,-1]))
          rownames(lr.result) = NULL
          
          message("====> LRT results:")
          print(lr.result)
          scm$lrt = lr.result
        }
        
        
      } else {
        if(cor){
          if (!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ str(~ vm(', colnames(control)[2],
                                                ',K) + grp(g1), ~corh(2):vm(',
                                                colnames(control)[2],',K)) + str(~ vm(',
                                                colnames(control)[2],',K):',
                                                colnames(control)[5],'+grp(g1):',
                                                colnames(control)[5], ', ~diag(2):id(',
                                                colnames(control)[5], '):vm(',
                                                colnames(control)[2], ',K))'))
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
          }
          
        }else{
          if (!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K)+ grp(g1) + vm(',
                                                colnames(control)[2],',K):',
                                                colnames(control)[5],' + grp(g1):',
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
          
        }
        
        if (any(na.exclude(scm$vparameters.pc > 1))) scm = .up.mod(scm)
        
        if(lrtest & scm$converge){
          message("====> Starting likelihood ratio tests")
          
          if (!is.null(K)){
            bench = suppressWarnings({
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+ vm(', colnames(control)[2],
                                            ',K)+ grp(g1) + vm(',
                                            colnames(control)[2],',K):',
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
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K)+ grp(g1)','+grp(g1):',
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
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K)+ grp(g1) + vm(',
                                                colnames(control)[2],',K):',
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
                                         paste0('+ grp(g1)',' + grp(g1):',
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
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K) + vm(',
                                                colnames(control)[2],',K):',
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
            
          }else{
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
                                                '+ grp(g1)','+ grp(g1):',
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
                                         paste0('+ grp(g1) ','+ grp(g1):',
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
            
          }
          
          
          lr.result = cbind(effect = c("DGE", 'IGE', 'DGE:age', 'IGE:age'), 
                            rbind(lr.d[,-1], lr.i[,-1], lr.dint[,-1], lr.iint[,-1]))
          rownames(lr.result) = NULL
          
          message("====> LRT results:")
          print(lr.result)
          scm$lrt = lr.result
        }
        
      }
    }
    
    
    # Single age====
    
  }else{
    if(inherits(prep.out, "comprepfor") && control[,6] > 1) ### Multi-areas ------------
    {
      input <<- input[order(input[,names(control)[6]],
                            input[,names(control)[3]], 
                            input[,names(control)[4]]),]
      
      if(spatial){
        if(cor){
          if(!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ str(~vm(', colnames(control)[2],
                                                ',K)+ grp(g1), ~corh(2):vm(',
                                                colnames(control)[2],',K))'))
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
          }
          
        }else{
          if(!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K)+ grp(g1)'))
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
                                         paste0('+ ', colnames(control)[2],
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
          
          
        }
        if (any(na.exclude(scm$vparameters.pc > 1))) scm = .up.mod(scm)
        
        if(lrtest & scm$converge){
          message("====> Starting likelihood ratio tests")
          
          if(!is.null(K)){
            bench = suppressWarnings({
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+ vm(', colnames(control)[2],
                                            ', K)+ grp(g1)'))
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
                                     paste0('+vm(', colnames(control)[2], ",K)"))
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
          }else{
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
          }
          
          lr.result = cbind(effect = c("DGE", 'IGE'), 
                            rbind(lr.d[,-1],lr.i[,-1]))
          rownames(lr.result) = NULL
          scm$lrt = lr.result
          message("====> LRT results:")
          print(lr.result)
        }
        
        
      }else{
        if(cor){
          if(!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ str(~vm(', colnames(control)[2],
                                                ',K)+ grp(g1), ~corh(2):vm(',
                                                colnames(control)[2],',K))'))
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
          }
          
        }else{
          if(!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+vm(', colnames(control)[2],
                                                ',K)+ grp(g1)'))
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
          
          
          if (any(na.exclude(scm$vparameters.pc > 1))) scm = .up.mod(scm)
        }
        if(lrtest & scm$converge){
          message("====> Starting likelihood ratio tests")
          
          if(!is.null(K)){
            bench = suppressWarnings({
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+vm(', colnames(control)[2],
                                            ',K)+ grp(g1)'))
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
                                     paste0('+vm(', colnames(control)[2],",K)"))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~dsum(~id(units)|',
                               colnames(control)[6],')')
                             ),
                             data = input, ...)
            )
            })
            
          }else{
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
          }
          
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
      input <<- input[order(input[,names(control)[3]], 
                            input[,names(control)[4]]),]
      
      if(spatial){
        if(cor){
          
          if(!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ str(~vm(', colnames(control)[2],
                                                ',K)+ grp(g1), ~corh(2):vm(',
                                                colnames(control)[2],',K))'))
                                 ),
                                 group = list(g1 = 1:control[,2]),
                                 residual = stats::as.formula(paste0(
                                   '~ar1(', colnames(control)[3],'):ar1(',
                                   colnames(control)[4],')')
                                 ),
                                 data = input, ...)
          }else{
            scm = asreml::asreml(
              fixed = fixed,
              random = stats::as.formula(paste(
                paste(as.character(random), collapse = ""),
                paste0(
                  '+ str(~',
                  colnames(control)[2],
                  '+ grp(g1), ~corh(2):id(',
                  colnames(control)[2],
                  '))'
                )
              )),
              group = list(g1 = 1:control[, 2]),
              residual = stats::as.formula(paste0(
                '~ar1(', colnames(control)[3], '):ar1(', colnames(control)[4], ')'
              )),
              data = input,
              ...
            )
          }
          
        }else{
          if(!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+vm(', colnames(control)[2],
                                                ',K)+ grp(g1)'))
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
          
        }
        
        if (any(na.exclude(scm$vparameters.pc > 1))) scm = .up.mod(scm)
        
        if(lrtest & scm$converge){
          message("====> Starting likelihood ratio tests")
          
          if(!is.null(K)){
            bench = suppressWarnings({
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+ vm(', colnames(control)[2],
                                            ',K)+ grp(g1)'))
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
                                     paste0('+ vm(', colnames(control)[2], ",K)"))
                             ),
                             group = list(g1 = 1:control[,2]),
                             residual = stats::as.formula(paste0(
                               '~ar1(', colnames(control)[3],'):ar1(',
                               colnames(control)[4],')')
                             ),
                             data = input, ...)
            )
            })
            
            
          }else{
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
            
          }
          
          
          lr.result = cbind(effect = c("DGE", 'IGE'), 
                            rbind(lr.d[,-1],lr.i[,-1]))
          rownames(lr.result) = NULL
          scm$lrt = lr.result
          message("====> LRT results:")
          print(lr.result)
        }
        
        
      }else{
        
        if(cor){
          if(!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ str(~ vm(', colnames(control)[2],
                                                ',K)+ grp(g1), ~corh(2):vm(',
                                                colnames(control)[2],',K))'))
                                 ),
                                 group = list(g1 = 1:control[,2]),
                                 data = input, ...)
          }else{
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ str(~', colnames(control)[2],
                                                '+ grp(g1), ~corh(2):id(',
                                                colnames(control)[2],'))'))
                                 ),
                                 group = list(g1 = 1:control[,2]),
                                 data = input, ...)
          }
          
        }else{
          if(!is.null(K)){
            scm = asreml::asreml(fixed = fixed,
                                 random = stats::as.formula(
                                   paste(paste(as.character(random), collapse = ""), 
                                         paste0('+ vm(', colnames(control)[2],
                                                ',K)+ grp(g1)'))
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
          
        }
        
        if (any(na.exclude(scm$vparameters.pc > 1))) scm = .up.mod(scm)
        
        if(lrtest & scm$converge){
          message("====> Starting likelihood ratio tests")
          
          if(!is.null(K)){
            bench = suppressWarnings({
              asreml::asreml(fixed = fixed,
                             random = stats::as.formula(
                               paste(paste(as.character(random), collapse = ""), 
                                     paste0('+vm(', colnames(control)[2],
                                            ',K)+ grp(g1)'))
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
                                     paste0('+vm(', colnames(control)[2],",K)"))
                             ),
                             group = list(g1 = 1:control[,2]),
                             data = input, ...)
            )
            })
          }else{
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
          }
          
          
          lr.result = cbind(effect = c("DGE", 'IGE'), 
                            rbind(lr.d[,-1],lr.i[,-1]))
          rownames(lr.result) = NULL
          scm$lrt = lr.result
          message("====> LRT results:")
          print(lr.result)
        }
      }
    }
  }
  
  class(scm) = c("compmod", class(scm))
  if(multi.age) attr(scm,"multi.age") = TRUE
  if(!is.null(K)) attr(scm, "treats") = attr(K, 'rowNames') else attr(scm, "treats") = levels(input[,colnames(control)[2]])
  if(!is.null(K)) remove(control, fixed, random, input, K, envir = .GlobalEnv) else remove(control, fixed, random, input, envir = .GlobalEnv)
  
  return(scm)
}

