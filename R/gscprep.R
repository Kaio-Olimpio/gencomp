gscprep<- function(data, gen, block, area = 1, ind, trait, row, col, dist.row, 
                   dist.col, method, n.dec = 2){
  
  x = data
  colnames(x)[which(colnames(x) == gen)] = 'trat'
  colnames(x)[which(colnames(x) == row)] = 'row'
  colnames(x)[which(colnames(x) == col)] = 'col'
  colnames(x)[which(colnames(x) == trait)] = 'trait'
  colnames(x)[which(colnames(x) == block)] = 'rep'
  colnames(x)[which(colnames(x) == ind)] = 'plot'
  colnames(x)[which(colnames(x) == area)] = 'area'
  
  if(!is.factor(x$trat)) x$trat = as.factor(x$trat)
  if(!is.factor(x$rep)) x$rep = as.factor(x$rep)
  if(!is.factor(x$plot)) x$plot = as.factor(x$plot)
  if(!is.factor(x$area)) x$area = as.factor(x$area)
  if(!is.numeric(x$row)) x$row = as.numeric(x$row)
  if(!is.numeric(x$col)) x$col = as.numeric(x$col)
  
  x = x[order(x$area,x$row,x$col),]
  p = dist.col/dist.row
  dist.diag = sqrt(dist.row^2 + dist.col^2)
  
  Z = list()
  z <- matrix(0, nrow(x), nlevels(x$trat), dimnames = list(1:nrow(x), levels(x$trat))) 
  
  fd = NULL
  fc = NULL
  fr = NULL
  nc = NULL
  nr = NULL
  nd = NULL
  w = list()
  for (i in 1:(nrow(z))) {
    n_r = n_rd = 0
    n_c = n_cd = 0 
    n_d = n_dd = 0 
    gen.diag = x$trat[x$row %in% (x$row[i] + c(-1, +1)) & x$col %in% (x$col[i] + c(-1, +1))]
    rep.diag = x$rep[x$row %in% (x$row[i] + c(-1, +1)) & x$col %in% (x$col[i] + c(-1, +1))]
    plot.diag = x$plot[x$row %in% (x$row[i] + c(-1, +1)) & x$col %in% (x$col[i] + c(-1, +1))]
    area.diag = x$area[x$row %in% (x$row[i] + c(-1, +1)) & x$col %in% (x$col[i] + c(-1, +1))]
    trt_d = NULL
    for (k in 1:length(gen.diag)) {
      if (is.na(x[which(x$trat %in% gen.diag[k] & 
                        x$rep %in% rep.diag[k] & 
                        x$plot %in% plot.diag[k]), "trait"]) | 
          x[i,'area'] != area.diag[k]){
        
        trt_d[k] = NA
      }else{
        n_d = n_d + 1
        trt_d[k] = x[which(x$trat %in% gen.diag[k] & 
                             x$rep %in% rep.diag[k] & 
                             x$plot %in% plot.diag[k]), "trait"]
      }
    }
    gen.col = x$trat[x$col == x$col[i] & x$row %in% (x$row[i] + c(-1,+1))]
    rep.col = x$rep[x$col == x$col[i] & x$row %in% (x$row[i] + c(-1,+1))]
    plot.col = x$plot[x$col == x$col[i] & x$row %in% (x$row[i] + c(-1,+1))]
    area.col = x$area[x$col == x$col[i] & x$row %in% (x$row[i] + c(-1,+1))]
    trt_c = NULL
    for (k in 1:length(gen.col)) {
      if (is.na(x[which(x$trat %in% gen.col[k] & 
                        x$rep %in% rep.col[k] & 
                        x$plot %in% plot.col[k]), "trait"]) | 
          x[i,'area'] != area.col[k]){
        n_cd = n_cd + 1
        trt_c[k] = NA
      }else{
        n_c = n_c + 1
        trt_c[k] = x[which(x$trat %in% gen.col[k] & 
                             x$rep %in% rep.col[k] & 
                             x$plot %in% plot.col[k]), "trait"]
      }
    }
    gen.row = x$trat[x$row == x$row[i] & x$col %in% (x$col[i] + c(-1, +1))]
    rep.row = x$rep[x$row == x$row[i] & x$col %in% (x$col[i] + c(-1, +1))]
    plot.row = x$plot[x$row == x$row[i] & x$col %in% (x$col[i] + c(-1, +1))]
    area.row = x$area[x$row == x$row[i] & x$col %in% (x$col[i] + c(-1, +1))]
    trt_r = NULL
    for (k in 1:length(gen.row)) {
      if (is.na(x[which(x$trat %in% gen.row[k] & 
                        x$rep %in% rep.row[k] & 
                        x$plot %in% plot.row[k]), "trait"]) | 
          x[i,'area'] != area.row[k]){
        n_rd = n_rd + 1
        trt_r[k] = NA
      }else{
        n_r = n_r + 1
        trt_r[k] = x[which(x$trat %in% gen.row[k] & 
                             x$rep %in% rep.row[k] & 
                             x$plot %in% plot.row[k]), "trait"]
      }
    }
    
    w[[i]] = data.frame(
      trat = x[which(x$col == x$col[i] & x$row == x$row[[i]]),'trat'],
      row = x$row[i],
      col = x$col[i],
      y_focal = x[which(x$col == x$col[i] & x$row == x$row[[i]]),'trait'],
      y_row = mean(trt_r, na.rm = T),
      #y_row_dist = mean(trt_r/dist.row, na.rm = T),
      n_row = n_r,
      y_col = mean(trt_c, na.rm = T),
      #y_col_dist = mean(trt_c/dist.col, na.rm = T),
      n_col = n_c,
      y_diag = mean(trt_d, na.rm = T),
      #y_diag_dist = mean(trt_d/dist.diag, na.rm = T),
      n_diag = n_d,
      y_neigh = mean(c(trt_r, trt_c, trt_d), na.rm = T) 
      #y_neigh_dist = mean(c(trt_r/dist.row, trt_c/dist.col, trt_d/dist.diag), na.rm = T)
    )
    
    if(method == "SK"){
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
    } else if(method == "CC"){
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
    }else if(method == "MU"){
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
  }
  
  # CIF
  cif = mean(fr, na.rm=T)*mean(nr) + mean(fc,na.rm=T)*mean(nc) + 
    mean(fd,na.rm=T)*mean(nd)
  
  #neigh_check
  w = do.call(rbind, w)
  
  #Field
  # x$trat_cod = ifelse(is.na(x$trait), paste0(x$trat,'*'), paste(x$trat))
  # field = stats::reshape(data = x[,c('row', 'col', 'trat_cod')], direction = 'wide', 
  #                        v.names = 'trat_cod', timevar = 'col', idvar = 'row')
  # rownames(field) = field$row; field = field[,-1]
  # colnames(field) = sub('trat_cod.', '', colnames(field))
  
  #Input
  data = data[order(data[,area], data[,row], data[,col]),]
  
  input = data.frame(cbind(z, data))
  input = input[order(input[,area], input[,row], input[,col]),]
  
  if(!is.factor(input[, gen])) input[, gen]= as.factor(input[, gen])
  if(!is.factor(input[, block])) input[, block] = as.factor(input[, block])
  if(!is.factor(input[, area])) input[, area] = as.factor(input[, area])
  if(!is.factor(input[, row])) input[, row] = as.factor(input[, row])
  if(!is.factor(input[, col])) input[, col] = as.factor(input[, col])
  
  Z = list(Z = z, CIF = cif, neigh_check = w, data = input)
  
  return(Z)
}