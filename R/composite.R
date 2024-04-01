##' @title Simulates a clonal composite
##'
##' @description
##' This function simulates clonal composites using outputs of a genetic-spatial
##' competition model fitted using [competition::asr()]
##' 
##' @param prep.out A `comprep` object.
##' @param model An `asreml` object, preferably obtained using the [competition::asr()] function.
##' @param resp.out A `comresp` object.
##' @param d.row.col A vector of size two. The first element contain the distance between
##' rows, and second the distance between columns of the simulated grid.
##' @param d.weight A logical value. If `TRUE` (default) the predicted mean 
##' of each plant in the grid will be weighted by the inverse of the distance between rows,
##' columns and diagonals.
##' @param selected A vector with the names of the clones selected to compose the clonal composite. 
##' The names must be present in the `comresp` object.
##' @param nsim An integer defining the number grid simulations. If `nsim > 1`, 
##' the function will estimate the 95% confidence interval of the predicted means using 
##' a bootstrap process. Defaults to 10.
##' @param verbose A logical value. If `TRUE`, shows the function progress. Defaults to `FALSE`.
##' 
##' @returns The function returns a data.frame with the predicted mean of each clone 
##' and their respective 95% confidence interval. The composite performance is obtained by
##' averaging the predicted mean of all clones.
##' 
##' @details
##' Considering the direct (DGE) and indirect genetic effects (IGE) of the selected clones, 
##' the function simulates grids. Clones are positioned differently in each simulation, 
##' which enables the modification of focal tree-neighbour dynamics. In each simulation, 
##' the expected mean of each clone is predicted using the following equation \insertCite{ferreira_novel_2023}{competition}:
##' \deqn{\hat{y}_{ij} = \mu + d_i + \sum^n_{i \neq j}{c_j}}
##' where \eqn{d_i} is the DGE of the i<sup>th</sup> focal tree, and 
##' \eqn{c_j} is the IGE of the j<sup>th</sup> neighbour. If `d.weight = TRUE`, the
##' IGE is divided by the distance between the focal tree and its neighbours:
##' \deqn{\hat{y}_{ij} = \mu + d_i + \sum^n_{i \neq j}{\frac{1}{dist_{ij}} \times c_j}}
##' 
##' @references 
##' \insertAllCited{}
##'
##' @seealso  [competition::prep], [competition::asr], [competition::resp]
##' 
##' @importFrom Rdpack reprompt
##' @importFrom stats quantile model.matrix model.matrix.lm
##' 
##' @export
##' 
##' @examples
##' \donttest{
##'  comp_mat = prep(data = euca, gen = 'clone', repl = 'block', area = 'area', 
##'                       ind = 'tree', age = 'age', row = 'row', col = 'col', 
##'                       dist.col = 3, dist.row = 2, trait = 'mai', method = 'SK',
##'                       n.dec = 3, verbose = TRUE)
##'  
##'  model = asr(prep.out = comp_mat, 
##'                   fixed = mai~ age, 
##'                   random = ~ block:age, 
##'                   cor = TRUE, maxit = 50)
##'                   
##'  results = resp(prep.out = comp_mat, model = model, weight.tgv = FALSE)
##'  
##'  cc = composite(prep.out = comp_mat, model = model, resp.out = results,
##'                 d.row.col = c(3,3), d.weight = TRUE, nsim = 10, verbose = TRUE, 
##'                 selected = results$blups$main[order(results$blups$main$TGV, 
##'                                                     decreasing = TRUE),1][1:10])
##'  }


composite = function(prep.out, model, resp.out, d.row.col, d.weight = TRUE, 
                     selected, nsim = 10, verbose = TRUE) {
  
  ## Objects from previous functions
  prep.out = prep.out
  resp.out = resp.out
  model = model
  
  # Distances
  drow = d.row.col[1]
  dcol = d.row.col[2]
  ddiag = sqrt(drow^2 + dcol^2)
  
  # New dataset: main effects
  dat = resp.out$blups$main
  dat = droplevels(dat[which(dat[,1] %in% selected),]); rownames(dat) = NULL
  ngen = length(unique(dat[,1]))
  mu = model$coefficients$fixed[grep('Intercept', rownames(model$coefficients$fixed))]
  
  # Simulating the grid: 1000 rows and 1000 columns
  grid = list()
  for (i in 1:nsim) {
    simu = data.frame(
      gen = sample(x = dat[,1], size = 100*100, replace = TRUE),
      row1 = sample(x = dat[,1], size = 100*100, replace = TRUE),
      row2 = sample(x = dat[,1], size = 100*100, replace = TRUE),
      col1 = sample(x = dat[,1], size = 100*100, replace = TRUE),
      col2 = sample(x = dat[,1], size = 100*100, replace = TRUE),
      diag1 = sample(x = dat[,1], size = 100*100, replace = TRUE),
      diag2 = sample(x = dat[,1], size = 100*100, replace = TRUE), 
      diag3 = sample(x = dat[,1], size = 100*100, replace = TRUE), 
      diag4 = sample(x = dat[,1], size = 100*100, replace = TRUE)
    )
    
    # Edges of the rectangle: 4
    simu[sample(nrow(simu), 4), c('row2', 'col2', 'diag2', 'diag3', 'diag4')] = NA
    
    # Horizontal sides of the rectangle: 1000 x 2
    simu[sample(nrow(simu), 200), c('row2', 'diag3', 'diag4')] = NA
    
    # Vertical sides of the rectangle: 1000 x 2
    simu[sample(nrow(simu), 200), c('col2', 'diag3', 'diag4')] = NA
    
    grid[[i]] = simu
  }
  
  if(verbose) cat("1. Simulated the grids", fill = TRUE)
  
  grid = lapply(grid, function(x){
    
    if(d.weight){
      dge = stats::model.matrix(~-1 + gen, data = x)
      colnames(dge) = sub('gen', '', colnames(dge))
      dge = dge %*% as.matrix(dat[match(dat[,1], colnames(dge)),'DGE'])
      
      ige.row1 = stats::model.matrix(~-1 + row1, data = x)
      colnames(ige.row1) = sub('row1', '', colnames(ige.row1))
      ige.row2 = stats::model.matrix.lm(~-1 + row2, data = x, na.action = 'na.pass')
      colnames(ige.row2) = sub('row2', '', colnames(ige.row2))
      ige.row2 = ifelse(is.na(ige.row2), 0, ige.row2)
      ige.row = ige.row1 + ige.row2; rm(ige.row1, ige.row2)
      ige.row = ige.row %*% (as.matrix(dat[match(dat[,1], colnames(ige.row)),'IGE'])/drow)
      
      ige.col1 = stats::model.matrix(~-1 + col1, data = x)
      colnames(ige.col1) = sub('col1', '', colnames(ige.col1))
      ige.col2 = stats::model.matrix.lm(~-1 + col2, data = x, na.action = 'na.pass')
      colnames(ige.col2) = sub('col2', '', colnames(ige.col2))
      ige.col2 = ifelse(is.na(ige.col2), 0, ige.col2)
      ige.col = ige.col1 + ige.col2; rm(ige.col1, ige.col2)
      ige.col = ige.col %*% (as.matrix(dat[match(dat[,1], colnames(ige.col)),'IGE'])/dcol)
      
      ige.diag1 = stats::model.matrix(~-1 + diag1, data = x)
      colnames(ige.diag1) = sub('diag1', '', colnames(ige.diag1))
      ige.diag2 = stats::model.matrix.lm(~-1 + diag2, data = x, na.action = 'na.pass')
      colnames(ige.diag2) = sub('diag2', '', colnames(ige.diag2))
      ige.diag2 = ifelse(is.na(ige.diag2), 0, ige.diag2)
      ige.diag3 = stats::model.matrix.lm(~-1 + diag3, data = x, na.action = 'na.pass')
      colnames(ige.diag3) = sub('diag3', '', colnames(ige.diag3))
      ige.diag3 = ifelse(is.na(ige.diag3), 0, ige.diag3)
      ige.diag4 = stats::model.matrix.lm(~-1 + diag4, data = x, na.action = 'na.pass')
      colnames(ige.diag4) = sub('diag4', '', colnames(ige.diag4))
      ige.diag4 = ifelse(is.na(ige.diag4), 0, ige.diag4)
      ige.diag = ige.diag1 + ige.diag2 + ige.diag3 + ige.diag4; rm(ige.diag1, ige.diag2, ige.diag3, ige.diag4)
      ige.diag = ige.diag %*% (as.matrix(dat[match(dat[,1], colnames(ige.diag)),'IGE'])/ddiag)
      
      data.frame(gen = x$gen, y.pred = mu + dge + ige.row + ige.col + ige.diag)
      
    }else{
      dge = stats::model.matrix(~-1 + gen, data = x)
      colnames(dge) = sub('gen', '', colnames(dge))
      dge = dge %*% as.matrix(dat[match(dat[,1], colnames(dge)),'DGE'])
      
      ige.row1 = stats::model.matrix(~-1 + row1, data = x)
      colnames(ige.row1) = sub('row1', '', colnames(ige.row1))
      ige.row2 = stats::model.matrix.lm(~-1 + row2, data = x, na.action = 'na.pass')
      colnames(ige.row2) = sub('row2', '', colnames(ige.row2))
      ige.row2 = ifelse(is.na(ige.row2), 0, ige.row2)
      ige.row = ige.row1 + ige.row2; rm(ige.row1, ige.row2)
      ige.row = ige.row %*% as.matrix(dat[match(dat[,1], colnames(ige.row)),'IGE'])
      
      ige.col1 = stats::model.matrix(~-1 + col1, data = x)
      colnames(ige.col1) = sub('col1', '', colnames(ige.col1))
      ige.col2 = stats::model.matrix.lm(~-1 + col2, data = x, na.action = 'na.pass')
      colnames(ige.col2) = sub('col2', '', colnames(ige.col2))
      ige.col2 = ifelse(is.na(ige.col2), 0, ige.col2)
      ige.col = ige.col1 + ige.col2; rm(ige.col1, ige.col2)
      ige.col = ige.col %*% as.matrix(dat[match(dat[,1], colnames(ige.col)),'IGE'])
      
      ige.diag1 = stats::model.matrix(~-1 + diag1, data = x)
      colnames(ige.diag1) = sub('diag1', '', colnames(ige.diag1))
      ige.diag2 = stats::model.matrix.lm(~-1 + diag2, data = x, na.action = 'na.pass')
      colnames(ige.diag2) = sub('diag2', '', colnames(ige.diag2))
      ige.diag2 = ifelse(is.na(ige.diag2), 0, ige.diag2)
      ige.diag3 = stats::model.matrix.lm(~-1 + diag3, data = x, na.action = 'na.pass')
      colnames(ige.diag3) = sub('diag3', '', colnames(ige.diag3))
      ige.diag3 = ifelse(is.na(ige.diag3), 0, ige.diag3)
      ige.diag4 = stats::model.matrix.lm(~-1 + diag4, data = x, na.action = 'na.pass')
      colnames(ige.diag4) = sub('diag4', '', colnames(ige.diag4))
      ige.diag4 = ifelse(is.na(ige.diag4), 0, ige.diag4)
      ige.diag = ige.diag1 + ige.diag2 + ige.diag3 + ige.diag4; rm(ige.diag1, ige.diag2, ige.diag3, ige.diag4)
      ige.diag = ige.diag %*% as.matrix(dat[match(dat[,1], colnames(ige.diag)),'IGE'])
      
      data.frame(gen = x$gen, y.pred = mu + dge + ige.row + ige.col + ige.diag)
    }
    
  })
  
  if(verbose) cat("2. Computed the total genotypic values", fill = TRUE)
  
  
  predgen = do.call(rbind, apply(do.call(cbind, lapply(grid, function(x) tapply(x$y.pred, x$gen, mean))),
        1, function(y){
          
          ## Bootstrap for computing the mean and the 95% confidence interval
          
          temp = vector()
          niter = 10000
          i = 0
          repeat{
            temp = c(temp, mean(sample(y, replace = TRUE)))
            i = i+1
            if(i == niter) break
          }
          
          CI = stats::quantile(sort(temp), prob = c(0.05, 0.95))
          data.frame(y.pred = mean(temp), CI_0.05 = CI[1], CI_0.95 = CI[2], 
                     row.names = NULL)
        }))
  
  if(verbose) cat("3. Estimated the predicted mean and its confidence intervals", fill = TRUE)
  
  # predgen$gen = rownames(predgen);rownames(predgen) = NULL;predgen = predgen[,c(4,1,2,3)]
  
  
  message(paste("\n The means were predicted considering an area of", 
                (dcol * drow * 100 * 100)/10000), ' ha')
  
  return(predgen)

}



