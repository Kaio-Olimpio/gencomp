composite = function(gscprep.out, model, gscresp.out, d.row.col, d.weight = TRUE, 
                     selected, nsim = 10, verbose = TRUE) {
  
  ## Objects from previous functions
  gscprep.out = gscprep.out
  gscresp.out = gscresp.out
  model = model
  
  # Distances
  drow = d.row.col[1]
  dcol = d.row.col[2]
  ddiag = sqrt(drow^2 + dcol^2)
  
  # New dataset: main effects
  dat = gscresp.out$blups$main
  dat = dat[which(dat[,1] %in% selected),]
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
      dge = model.matrix(~-1 + gen, data = x)
      colnames(dge) = sub('gen', '', colnames(dge))
      dge = dge %*% as.matrix(dat[match(dat[,1], colnames(dge)),'DGE'])
      
      ige.row1 = model.matrix(~-1 + row1, data = x)
      colnames(ige.row1) = sub('row1', '', colnames(ige.row1))
      ige.row2 = model.matrix.lm(~-1 + row2, data = x, na.action = 'na.pass')
      colnames(ige.row2) = sub('row2', '', colnames(ige.row2))
      ige.row2 = ifelse(is.na(ige.row2), 0, ige.row2)
      ige.row = ige.row1 + ige.row2; rm(ige.row1, ige.row2)
      ige.row = ige.row %*% (as.matrix(dat[match(dat[,1], colnames(ige.row)),'IGE'])/drow)
      
      ige.col1 = model.matrix(~-1 + col1, data = x)
      colnames(ige.col1) = sub('col1', '', colnames(ige.col1))
      ige.col2 = model.matrix.lm(~-1 + col2, data = x, na.action = 'na.pass')
      colnames(ige.col2) = sub('col2', '', colnames(ige.col2))
      ige.col2 = ifelse(is.na(ige.col2), 0, ige.col2)
      ige.col = ige.col1 + ige.col2; rm(ige.col1, ige.col2)
      ige.col = ige.col %*% (as.matrix(dat[match(dat[,1], colnames(ige.col)),'IGE'])/dcol)
      
      ige.diag1 = model.matrix(~-1 + diag1, data = x)
      colnames(ige.diag1) = sub('diag1', '', colnames(ige.diag1))
      ige.diag2 = model.matrix.lm(~-1 + diag2, data = x, na.action = 'na.pass')
      colnames(ige.diag2) = sub('diag2', '', colnames(ige.diag2))
      ige.diag2 = ifelse(is.na(ige.diag2), 0, ige.diag2)
      ige.diag3 = model.matrix.lm(~-1 + diag3, data = x, na.action = 'na.pass')
      colnames(ige.diag3) = sub('diag3', '', colnames(ige.diag3))
      ige.diag3 = ifelse(is.na(ige.diag3), 0, ige.diag3)
      ige.diag4 = model.matrix.lm(~-1 + diag4, data = x, na.action = 'na.pass')
      colnames(ige.diag4) = sub('diag4', '', colnames(ige.diag4))
      ige.diag4 = ifelse(is.na(ige.diag4), 0, ige.diag4)
      ige.diag = ige.diag1 + ige.diag2 + ige.diag3 + ige.diag4; rm(ige.diag1, ige.diag2, ige.diag3, ige.diag4)
      ige.diag = ige.diag %*% (as.matrix(dat[match(dat[,1], colnames(ige.diag)),'IGE'])/ddiag)
      
      data.frame(gen = x$gen, y.pred = mu + dge + ige.row + ige.col + ige.diag)
      
    }else{
      dge = model.matrix(~-1 + gen, data = x)
      colnames(dge) = sub('gen', '', colnames(dge))
      dge = dge %*% as.matrix(dat[match(dat[,1], colnames(dge)),'DGE'])
      
      ige.row1 = model.matrix(~-1 + row1, data = x)
      colnames(ige.row1) = sub('row1', '', colnames(ige.row1))
      ige.row2 = model.matrix.lm(~-1 + row2, data = x, na.action = 'na.pass')
      colnames(ige.row2) = sub('row2', '', colnames(ige.row2))
      ige.row2 = ifelse(is.na(ige.row2), 0, ige.row2)
      ige.row = ige.row1 + ige.row2; rm(ige.row1, ige.row2)
      ige.row = ige.row %*% as.matrix(dat[match(dat[,1], colnames(ige.row)),'IGE'])
      
      ige.col1 = model.matrix(~-1 + col1, data = x)
      colnames(ige.col1) = sub('col1', '', colnames(ige.col1))
      ige.col2 = model.matrix.lm(~-1 + col2, data = x, na.action = 'na.pass')
      colnames(ige.col2) = sub('col2', '', colnames(ige.col2))
      ige.col2 = ifelse(is.na(ige.col2), 0, ige.col2)
      ige.col = ige.col1 + ige.col2; rm(ige.col1, ige.col2)
      ige.col = ige.col %*% as.matrix(dat[match(dat[,1], colnames(ige.col)),'IGE'])
      
      ige.diag1 = model.matrix(~-1 + diag1, data = x)
      colnames(ige.diag1) = sub('diag1', '', colnames(ige.diag1))
      ige.diag2 = model.matrix.lm(~-1 + diag2, data = x, na.action = 'na.pass')
      colnames(ige.diag2) = sub('diag2', '', colnames(ige.diag2))
      ige.diag2 = ifelse(is.na(ige.diag2), 0, ige.diag2)
      ige.diag3 = model.matrix.lm(~-1 + diag3, data = x, na.action = 'na.pass')
      colnames(ige.diag3) = sub('diag3', '', colnames(ige.diag3))
      ige.diag3 = ifelse(is.na(ige.diag3), 0, ige.diag3)
      ige.diag4 = model.matrix.lm(~-1 + diag4, data = x, na.action = 'na.pass')
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
          
          CI = stats::quantile(temp, prob = c(0.05, 0.95))
          data.frame(y.pred = mean(temp), CI_0.05 = CI[1], CI_0.95 = CI[2], 
                     row.names = NULL)
        }))
  
  if(verbose) cat("3. Estimated the predicted mean and its confidence intervals", fill = TRUE)
  
  # predgen$gen = rownames(predgen);rownames(predgen) = NULL;predgen = predgen[,c(4,1,2,3)]
  
  
  message(paste("\n The means were predicted considering an area of", 
                (dcol * drow * 100 * 100)/10000), ' ha')
  
  return(predgen)

}



