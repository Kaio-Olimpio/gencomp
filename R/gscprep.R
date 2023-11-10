##'
##'
##'
##'
##'
##'
##'
##'
##' 
##' @export

gscprep<- function(data, gen, repl, row, col, ind, trait, dist.row, dist.col, 
                   method = 'SK', area = NULL, age = NULL, n.dec = 2, 
                   verbose = FALSE){
  
  # Messages and warnings
  if(length(unique(data[,ind])) != 
     length(unique(data[,row])) * length(unique(data[,col]))){
    warning("The trial does not have a rectangular grid!")
  }
  
  # Entry -------------------------------------------------------------------
  p = dist.col/dist.row
  dist.diag = sqrt(dist.row^2 + dist.col^2)
  
  x = data.frame(
    trat = as.factor(data[, gen]),
    repl = as.factor(data[, repl]),
    ind = as.factor(data[, ind]), 
    row = as.numeric(data[, row]),
    col = as.numeric(data[, col]),
    trait = as.numeric(data[, trait])
  )
  
  control = data.frame(
    trait = length(x$trait), 
    ngen = nlevels(x$trat),
    nrepl = nlevels(x$repl),
    nrow = length(unique(x$row)),
    ncol = length(unique(x$col))
  )
  
  colnames(control) = c(trait, gen, repl, row, col)
  
  if(is.null(age)){
    control$age = 0
  }else{
    control[, age] = length(unique(data[, age]))
  }
  
  if(is.null(area)){
    control$area = 0
  }else{
    control[, area] = length(unique(data[, area]))
  }
  
  if(!is.null(area)){ # If the trial was laid out in different areas (talhões)
    x$area = as.factor(data[, area])  
  } else {
    x$area = 1
  }
  if(is.null(age)) # A single age ----------------------
  {
    if(is.null(area)){
      x = x[order(x$row, x$col),]
    }else{
      x = x[order(x$area, x$row, x$col),]
    }
    
    Z = list()
    z <- matrix(0, nrow(x), nlevels(x$trat), dimnames = list(1:nrow(x), levels(x$trat))) 
    fd = NULL
    fc = NULL
    fr = NULL
    nc = NULL
    nr = NULL
    nd = NULL
    w = list()
    
    for (i in 1:nrow(z)) ### Z competition matrix -----------
    {
      n_r = n_rd = 0
      n_c = n_cd = 0 
      n_d = n_dd = 0 
      
      ## Diagonal --------------
      gen.diag = x$trat[x$row %in% (x$row[i] + c(-1, +1)) & x$col %in% (x$col[i] + c(-1, +1))]
      repl.diag = x$repl[x$row %in% (x$row[i] + c(-1, +1)) & x$col %in% (x$col[i] + c(-1, +1))]
      ind.diag = x$ind[x$row %in% (x$row[i] + c(-1, +1)) & x$col %in% (x$col[i] + c(-1, +1))]
      area.diag = x$area[x$row %in% (x$row[i] + c(-1, +1)) & x$col %in% (x$col[i] + c(-1, +1))]
      trt_d = NULL
      for (k in 1:length(gen.diag)) {
        if (is.na(x[which(x$trat %in% gen.diag[k] & 
                          x$repl %in% repl.diag[k] & 
                          x$ind %in% ind.diag[k]), "trait"]) | 
            x[i,'area'] != area.diag[k]){
          
          trt_d[k] = NA
        }else{
          n_d = n_d + 1
          trt_d[k] = x[which(x$trat %in% gen.diag[k] & 
                               x$repl %in% repl.diag[k] & 
                               x$ind %in% ind.diag[k]), "trait"]
        }
      }
      
      ## Column ---------------
      gen.col = x$trat[x$col == x$col[i] & x$row %in% (x$row[i] + c(-1,+1))]
      repl.col = x$repl[x$col == x$col[i] & x$row %in% (x$row[i] + c(-1,+1))]
      ind.col = x$ind[x$col == x$col[i] & x$row %in% (x$row[i] + c(-1,+1))]
      area.col = x$area[x$col == x$col[i] & x$row %in% (x$row[i] + c(-1,+1))]
      trt_c = NULL
      for (k in 1:length(gen.col)) {
        if (is.na(x[which(x$trat %in% gen.col[k] & 
                          x$repl %in% repl.col[k] & 
                          x$ind %in% ind.col[k]), "trait"]) | 
            x[i,'area'] != area.col[k]){
          n_cd = n_cd + 1
          trt_c[k] = NA
        }else{
          n_c = n_c + 1
          trt_c[k] = x[which(x$trat %in% gen.col[k] & 
                               x$repl %in% repl.col[k] & 
                               x$ind %in% ind.col[k]), "trait"]
        }
      }
      
      ## Row ------------------
      gen.row = x$trat[x$row == x$row[i] & x$col %in% (x$col[i] + c(-1, +1))]
      repl.row = x$repl[x$row == x$row[i] & x$col %in% (x$col[i] + c(-1, +1))]
      ind.row = x$ind[x$row == x$row[i] & x$col %in% (x$col[i] + c(-1, +1))]
      area.row = x$area[x$row == x$row[i] & x$col %in% (x$col[i] + c(-1, +1))]
      trt_r = NULL
      for (k in 1:length(gen.row)) {
        if (is.na(x[which(x$trat %in% gen.row[k] & 
                          x$repl %in% repl.row[k] & 
                          x$ind %in% ind.row[k]), "trait"]) | 
            x[i,'area'] != area.row[k]){
          n_rd = n_rd + 1
          trt_r[k] = NA
        }else{
          n_r = n_r + 1
          trt_r[k] = x[which(x$trat %in% gen.row[k] & 
                               x$repl %in% repl.row[k] & 
                               x$ind %in% ind.row[k]), "trait"]
        }
      }
      
      ### Neighbourhood check -------------
      w[[i]] = data.frame(
        trat = x[which(x$col == x$col[i] & x$row == x$row[[i]]),'trat'],
        row = x$row[i],
        col = x$col[i],
        y_focal = x[which(x$col == x$col[i] & x$row == x$row[[i]]),'trait'],
        y_row = mean(trt_r, na.rm = T),
        n_row = n_r,
        y_col = mean(trt_c, na.rm = T),
        n_col = n_c,
        y_diag = mean(trt_d, na.rm = T),
        n_diag = n_d,
        y_neigh = mean(c(trt_r, trt_c, trt_d), na.rm = T) 
      )
      
      ### Competition intensity methods -------------
      
      if(method == "SK"){ # Costa e Silva and Kerr
        if(n_c == 0 & n_r == 0 & n_d == 0){
          f_D = 0
          z[i, gen.diag[which(area.diag == x[i,'area'])]] = f_D
          z[i, gen.col[which(area.col == x[i,'area'])]] = f_C = round((f_D * sqrt(1 + p ^ 2)) / p, n.dec) #fc
          z[i, gen.row[which(area.row == x[i,'area'])]] = f_R = round(f_D * sqrt(1 + p ^ 2), n.dec)#fr
          
          fd[i] = f_D
          fc[i] = f_C
          fr[i] = f_R
          nc[i] = n_c
          nr[i] = n_r
          nd[i] = n_d
        }else{
          f_D = round(p / (sqrt((n_r * p ^ 4) + (n_r * p ^ 2)
                                + (n_c * p ^ 2) + (n_d * p ^ 2) + n_c)), n.dec) 
          z[i, gen.diag[which(area.diag == x[i,'area'])]] = f_D
          z[i, gen.col[which(area.col == x[i,'area'])]] = f_C = round((f_D * sqrt(1 + p ^ 2)) / p, n.dec) #fc
          z[i, gen.row[which(area.row == x[i,'area'])]] = f_R = round(f_D * sqrt(1 + p ^ 2), n.dec)#fr
          
          fd[i] = f_D
          fc[i] = f_C
          fr[i] = f_R
          nc[i] = n_c
          nr[i] = n_r
          nd[i] = n_d
        }
      } else if(method == "CC"){  # Cappa and Cantet
        if(n_c == 0 & n_r == 0 & n_d == 0){
          z[i, gen.diag[which(area.diag == x[i,'area'])]] = f_D =  0
          z[i, gen.col[which(area.col == x[i,'area'])]] = f_C = 0  
          z[i, gen.row[which(area.row == x[i,'area'])]] = f_R = 0
          
          fd[i] = f_D
          fc[i] = f_C
          fr[i] = f_R
          nc[i] = n_c
          nr[i] = n_r
          nd[i] = n_d
        }else{
          z[i, gen.diag[which(area.diag == x[i,'area'])]] = f_D = round(1/sqrt(2*(n_r+n_c) + n_d), n.dec)
          z[i, gen.col[which(area.col == x[i,'area'])]] = f_C = round(sqrt(2/(2*(n_r+n_c) + n_d)), n.dec)
          z[i, gen.row[which(area.row == x[i,'area'])]] = f_R = round(sqrt(2/(2*(n_r+n_c) + n_d)), n.dec)
          
          fd[i] = f_D
          fc[i] = f_C
          fr[i] = f_R
          nc[i] = n_c
          nr[i] = n_r
          nd[i] = n_d
        }
      }else if(method == "MU"){ # Muir
        z[i, gen.diag[which(area.diag == x[i,'area'])]] = f_D = round((1/dist.diag), n.dec)
        z[i, gen.col[which(area.col == x[i,'area'])]] = f_C = round((1/dist.col), n.dec)
        z[i, gen.row[which(area.row == x[i,'area'])]] = f_R = round((1/dist.row), n.dec)
        
        fd[i] = f_D
        fc[i] = f_C
        fr[i] = f_R
        nc[i] = n_c
        nr[i] = n_r
        nd[i] = n_d
      }
      
      ### Progress bar
      if(verbose){
        cat('\r Running through the grid -->',
            sprintf(
              '%s%s|% 3s%%',
              strrep('=', round(i/nrow(z) * (options()$width - nchar('||100%')))),
              strrep(' ', options()$width -
                       round(i/nrow(z) * (options()$width - nchar('||100%'))) -
                       nchar('||100%')),
              round(i/nrow(z)*100)
            )
        )
      }
      
    }
    
    cif = mean(fr, na.rm=T)*mean(nr) + mean(fc,na.rm=T)*mean(nc) + 
      mean(fd,na.rm=T)*mean(nd)
    w = do.call(rbind, w)
    
    ### Entry for the model function
    
    if(is.null(area)){
      data = data[order(data[, row], data[, col]),]
    } else{
      data = data[order(data[, area], data[, row], data[, col]),]
    }
    
    input = data.frame(cbind(z, data))
    
    if(is.null(area)){
      input = input[order(input[, row], input[, col]),]
    } else{
      input = input[order(input[, area], input[, row], input[, col]),]
    }
    
    if(!is.factor(input[, gen])) input[, gen]= as.factor(input[, gen])
    if(!is.factor(input[, repl])) input[, repl] = as.factor(input[, repl])
    if(!is.factor(input[, row])) input[, row] = as.factor(input[, row])
    if(!is.factor(input[, col])) input[, col] = as.factor(input[, col])
    if(!is.null(area)){
      if(!is.factor(input[, area])) input[, area] = as.factor(input[, area])
    }
    
    Z = list(Z = z, CIF = cif, neigh_check = w, data = input, control = control)
    
    return(Z)
    
  } else # Several ages -----------------------
  { 
    x$age = as.factor(data[, age])  
    x = split(x, x$age)
    x = lapply(x, function(q){
      if(is.null(area)){
        q[order(q$row, q$col),]
      }else{
        q[order(q$area, q$row, q$col),]
      }
    })
    
    Z = lapply(x, function(q){
      z <- matrix(0, nrow(q), nlevels(q$trat), dimnames = list(1:nrow(q), levels(q$trat))) 
      fd = NULL
      fc = NULL
      fr = NULL
      nc = NULL
      nr = NULL
      nd = NULL
      w = list()
      
      for (i in 1:nrow(z)) ### Z competition matrix -----------
      {
        n_r = n_rd = 0
        n_c = n_cd = 0 
        n_d = n_dd = 0 
        
        ## Diagonal --------------
        gen.diag = q$trat[q$row %in% (q$row[i] + c(-1, +1)) & q$col %in% (q$col[i] + c(-1, +1))]
        repl.diag = q$repl[q$row %in% (q$row[i] + c(-1, +1)) & q$col %in% (q$col[i] + c(-1, +1))]
        ind.diag = q$ind[q$row %in% (q$row[i] + c(-1, +1)) & q$col %in% (q$col[i] + c(-1, +1))]
        area.diag = q$area[q$row %in% (q$row[i] + c(-1, +1)) & q$col %in% (q$col[i] + c(-1, +1))]
        trt_d = NULL
        for (k in 1:length(gen.diag)) {
          if (is.na(q[which(q$trat %in% gen.diag[k] & 
                            q$repl %in% repl.diag[k] & 
                            q$ind %in% ind.diag[k]), "trait"]) | 
              q[i,'area'] != area.diag[k]){
            
            trt_d[k] = NA
          }else{
            n_d = n_d + 1
            trt_d[k] = q[which(q$trat %in% gen.diag[k] & 
                                 q$repl %in% repl.diag[k] & 
                                 q$ind %in% ind.diag[k]), "trait"]
          }
        }
        
        ## Column ---------------
        gen.col = q$trat[q$col == q$col[i] & q$row %in% (q$row[i] + c(-1,+1))]
        repl.col = q$repl[q$col == q$col[i] & q$row %in% (q$row[i] + c(-1,+1))]
        ind.col = q$ind[q$col == q$col[i] & q$row %in% (q$row[i] + c(-1,+1))]
        area.col = q$area[q$col == q$col[i] & q$row %in% (q$row[i] + c(-1,+1))]
        trt_c = NULL
        for (k in 1:length(gen.col)) {
          if (is.na(q[which(q$trat %in% gen.col[k] & 
                            q$repl %in% repl.col[k] & 
                            q$ind %in% ind.col[k]), "trait"]) | 
              q[i,'area'] != area.col[k]){
            n_cd = n_cd + 1
            trt_c[k] = NA
          }else{
            n_c = n_c + 1
            trt_c[k] = q[which(q$trat %in% gen.col[k] & 
                                 q$repl %in% repl.col[k] & 
                                 q$ind %in% ind.col[k]), "trait"]
          }
        }
        
        ## Row ------------------
        gen.row = q$trat[q$row == q$row[i] & q$col %in% (q$col[i] + c(-1, +1))]
        repl.row = q$repl[q$row == q$row[i] & q$col %in% (q$col[i] + c(-1, +1))]
        ind.row = q$ind[q$row == q$row[i] & q$col %in% (q$col[i] + c(-1, +1))]
        area.row = q$area[q$row == q$row[i] & q$col %in% (q$col[i] + c(-1, +1))]
        trt_r = NULL
        for (k in 1:length(gen.row)) {
          if (is.na(q[which(q$trat %in% gen.row[k] & 
                            q$repl %in% repl.row[k] & 
                            q$ind %in% ind.row[k]), "trait"]) | 
              q[i,'area'] != area.row[k]){
            n_rd = n_rd + 1
            trt_r[k] = NA
          }else{
            n_r = n_r + 1
            trt_r[k] = q[which(q$trat %in% gen.row[k] & 
                                 q$repl %in% repl.row[k] & 
                                 q$ind %in% ind.row[k]), "trait"]
          }
        }
        
        ### Neighbourhood check -------------
        w[[i]] = data.frame(
          trat = q[which(q$col == q$col[i] & q$row == q$row[[i]]),'trat'],
          row = q$row[i],
          col = q$col[i],
          y_focal = q[which(q$col == q$col[i] & q$row == q$row[[i]]),'trait'],
          y_row = mean(trt_r, na.rm = T),
          n_row = n_r,
          y_col = mean(trt_c, na.rm = T),
          n_col = n_c,
          y_diag = mean(trt_d, na.rm = T),
          n_diag = n_d,
          y_neigh = mean(c(trt_r, trt_c, trt_d), na.rm = T) 
        )
        
        ### Competition intensity methods -------------
        
        if(method == "SK"){ # Costa e Silva and Kerr
          if(n_c == 0 & n_r == 0 & n_d == 0){
            f_D = 0
            z[i, gen.diag[which(area.diag == q[i,'area'])]] = f_D
            z[i, gen.col[which(area.col == q[i,'area'])]] = f_C = round((f_D * sqrt(1 + p ^ 2)) / p, n.dec) #fc
            z[i, gen.row[which(area.row == q[i,'area'])]] = f_R = round(f_D * sqrt(1 + p ^ 2), n.dec)#fr
            
            fd[i] = f_D
            fc[i] = f_C
            fr[i] = f_R
            nc[i] = n_c
            nr[i] = n_r
            nd[i] = n_d
          }else{
            f_D = round(p / (sqrt((n_r * p ^ 4) + (n_r * p ^ 2)
                                  + (n_c * p ^ 2) + (n_d * p ^ 2) + n_c)), n.dec) 
            z[i, gen.diag[which(area.diag == q[i,'area'])]] = f_D
            z[i, gen.col[which(area.col == q[i,'area'])]] = f_C = round((f_D * sqrt(1 + p ^ 2)) / p, n.dec) #fc
            z[i, gen.row[which(area.row == q[i,'area'])]] = f_R = round(f_D * sqrt(1 + p ^ 2), n.dec)#fr
            
            fd[i] = f_D
            fc[i] = f_C
            fr[i] = f_R
            nc[i] = n_c
            nr[i] = n_r
            nd[i] = n_d
          }
        } else if(method == "CC"){  # Cappa and Cantet
          if(n_c == 0 & n_r == 0 & n_d == 0){
            z[i, gen.diag[which(area.diag == q[i,'area'])]] = f_D =  0
            z[i, gen.col[which(area.col == q[i,'area'])]] = f_C = 0  
            z[i, gen.row[which(area.row == q[i,'area'])]] = f_R = 0
            
            fd[i] = f_D
            fc[i] = f_C
            fr[i] = f_R
            nc[i] = n_c
            nr[i] = n_r
            nd[i] = n_d
          }else{
            z[i, gen.diag[which(area.diag == q[i,'area'])]] = f_D = round(1/sqrt(2*(n_r+n_c) + n_d), n.dec)
            z[i, gen.col[which(area.col == q[i,'area'])]] = f_C = round(sqrt(2/(2*(n_r+n_c) + n_d)), n.dec)
            z[i, gen.row[which(area.row == q[i,'area'])]] = f_R = round(sqrt(2/(2*(n_r+n_c) + n_d)), n.dec)
            
            fd[i] = f_D
            fc[i] = f_C
            fr[i] = f_R
            nc[i] = n_c
            nr[i] = n_r
            nd[i] = n_d
          }
        }else if(method == "MU"){ # Muir
          z[i, gen.diag[which(area.diag == q[i,'area'])]] = f_D = round((1/dist.diag), n.dec)
          z[i, gen.col[which(area.col == q[i,'area'])]] = f_C = round((1/dist.col), n.dec)
          z[i, gen.row[which(area.row == q[i,'area'])]] = f_R = round((1/dist.row), n.dec)
          
          fd[i] = f_D
          fc[i] = f_C
          fr[i] = f_R
          nc[i] = n_c
          nr[i] = n_r
          nd[i] = n_d
        }
        
        ### Progress bar
        if(verbose){
          cat(paste0('\r Running through the grid (Age ', unique(q$age),') -->'),
              sprintf(
                '%s%s|% 3s%%',
                strrep('=', round(i/nrow(z) * (options()$width - nchar('||100%')))),
                strrep(' ', options()$width -
                         round(i/nrow(z) * (options()$width - nchar('||100%'))) -
                         nchar('||100%')),
                round(i/nrow(z)*100)
              )
          )
        }
      }
      
      cif = mean(fr, na.rm=T)*mean(nr) + mean(fc,na.rm=T)*mean(nc) + 
        mean(fd,na.rm=T)*mean(nd)
      w = do.call(rbind, w)
      
      ### Entry for the model function
      
      dat = data[data[, age] == q$age, ]
      
      if(is.null(area)){
        dat = dat[order(dat[, row], dat[, col]),]
      } else{
        dat = dat[order(dat[, area], dat[, row], dat[, col]),]
      }
      
      input = data.frame(cbind(z, dat))
      
      if(!is.factor(input[, age])) input[, age]= as.factor(input[, age])
      if(!is.factor(input[, gen])) input[, gen]= as.factor(input[, gen])
      if(!is.factor(input[, repl])) input[, repl] = as.factor(input[, repl])
      if(!is.factor(input[, row])) input[, row] = as.factor(input[, row])
      if(!is.factor(input[, col])) input[, col] = as.factor(input[, col])
      if(!is.null(area)){
        if(!is.factor(input[, area])) input[, area] = as.factor(input[, area])
      }
      
      list(Z = z, CIF = cif, neigh_check = w, data = input)
    })
    
    names(Z) = paste0("Age_", names(Z))
    Z$data = rbind(Z$Age_3$data, Z$Age_6$data)
    Z$control = control
    
    return(Z)
  }
}

#Field
# x$trat_cod = ifelse(is.na(x$trait), paste0(x$trat,'*'), paste(x$trat))
# field = stats::reshape(data = x[,c('row', 'col', 'trat_cod')], direction = 'wide', 
#                        v.names = 'trat_cod', timevar = 'col', idvar = 'row')
# rownames(field) = field$row; field = field[,-1]
# colnames(field) = sub('trat_cod.', '', colnames(field))