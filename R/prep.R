##' @title Preparations to fit a genetic-spatial competition model
##' 
##' @description
##' This function builds the genetic competition matrix (\eqn{\mathbf{Z}_c}), prepare the 
##' dataset to be used to fit the model and provide summaries of the trial. It also
##' computes the competition intensity factor.
##' 
##' @param data A data frame containing the phenotypic data.
##' @param gen,repl,row,col,ind,trait A string. The name of the columns that correspond
##' to the genotype, replicate, row, column, individual and trait information,
##' respectively. `ind` can also distinguish the plots, if there are multiple plants per plot.
##' @param dist.row,dist.col An integer. The distance between rows and columns in the trial.
##' @param area A string. The name of the column that corresponds to the area information. 
##' Valid if the trial has non-contiguous blocks, for e.g., blocks 1 and 2 in area 1, and
##' blocks 3 and 4 in area 2. `NULL` (default) otherwise.
##' @param age A string. The name of the column that corresponds to the age information.
##' Necessary for fitting a multi-age model using [competition::asr()]. `NULL` (default)
##' otherwise.
##' @param method A string. The method for computing the competition intensity in \eqn{\mathbf{Z}_c}. 
##' It has three options: "MU" for the method proposed by \insertCite{muir_incorporation_2005;textual}{competition}, 
##' "CC" for the method proposed by \insertCite{cappa_direct_2008;textual}{competition}, and."SK" for the 
##' method proposed by \insertCite{costa_e_silva_accounting_2013;textual}{competition}.
##' See Details for more information on these methods.
##' @param n.dec An integer. The number of decimal digits to be displayed in \eqn{\mathbf{Z}_c}. 
##' Defaults to 2.
##' @param verbose A logical value. If `TRUE`, a progress bar will be displayed in the 
##' console. Defaults to `FALSE`.
##' 
##' 
##' @return The function returns:
##' \itemize{
##' \item \code{Z} : The \eqn{\mathbf{Z}_c} matrix
##' \item \code{data} : A data frame composed of the built \eqn{\mathbf{Z}_c} merged
##' with the dataset provided by the user.
##' \item \code{neigh_check} : A data frame containing the phenotypic records of each
##' focal plant and its neighbours.
##' \item \code{control} : A data frame containing the number of phenotypic records, 
##' genotypes, replicates, rows, columns, ages and areas in the dataset provided by the user.
##' }
##' If `age` is not `NULL`, the function provides the outputs described for each age.
##' 
##' @details
##' Three methods are available for estimating the competition intensity and building
##' the \eqn{\mathbf{Z}_c}, the genetic competition matrix: 
##' 
##' \itemize{\item \insertCite{muir_incorporation_2005;textual}{competition}: "MU"}
##' 
##' The average competition intensity is the inverse of the distance between the focal 
##' plant and its neighbours:
##' 
##' \deqn{f_{D_v} = \frac{1}{\sqrt{d_R^2 + d_C^2}}}
##' 
##' \deqn{f_{R_v} = \frac{1}{d_R}}
##' 
##' \deqn{f_{C_v} = \frac{1}{d_C}}
##'
##' where \eqn{f_{D_v}}, \eqn{f_{R_v}} and \eqn{f_{C_v}} are the average competition intensities 
##' in the diagonal, row, and column directions of the \eqn{v^{th}} clone, 
##' respectively; and \eqn{d_R} and \eqn{d_C} are the inter-row and inter-column distances. 
##'
##' \itemize{\item \insertCite{cappa_direct_2008;textual}{competition}: "CC"}
##' 
##' The average competition intensity depends on the number of neighbours in each direction
##' 
##' \deqn{f_{D_v} = \frac{1}{\sqrt{2(n_{C_v} + n_{R_v}) + n_{D_v}}}}
##' 
##' \deqn{f_{R_v} = \sqrt{\frac{2}{2(n_{C_v} + n_{R_v}) + n_{D_v}}}}
##' 
##' \deqn{f_{C_v} = \sqrt{\frac{2}{2(n_{C_v} + n_{R_v}) + n_{D_v}}}}
##' 
##' where \eqn{n_{D_v}}, \eqn{n_{R_v}} and \eqn{n_{C_v}} are the number of neighbours in the 
##' diagonal, row and column directions of the \eqn{v^{th}} clone, respectively.
##' Note that, in this case, it is assumed that the distance between rows and 
##' columns are the same.
##' 
##' \itemize{\item \insertCite{costa_e_silva_accounting_2013;textual}{competition}: "SK"}
##' 
##' The average competition intensity depends on both the distance between the focal 
##' tree and its neighbours, and the number of neighbours in each direction:
##'
##' \deqn{f_{D_v} = \frac{p}{\sqrt{(n_{R_v} p^4) + (n_{R_v} p^2) + (n_{C_v} p^2) + (n_{D_v} p^2) + n_{C_v}}}}
##' 
##' \deqn{f_{R_v} = f_{D_v} \sqrt{1 + p^2}}
##' 
##' \deqn{f_{C_v} = \frac{f_{D_v} \sqrt{1 + p^2}}{p}}
##' 
##' where \eqn{p = \frac{d_C}{d_R}}.
##' 
##' The overall competition intensity factor (\eqn{CIF}) is estimated by taking the mean of 
##' the competition intensities provided by the equations described above, and the
##' mean number of neighbours in each direction:
##'
##' \deqn{CIF = \overline{n_D} \overline{f_D} + \overline{n_R} \overline{f_R} + \overline{n_C} \overline{f_C}}
##' 
##' @references 
##' \insertAllCited{}
##' 
##' @importFrom Rdpack reprompt
##' 
##' @export
##' 
##' 
##' @examples
##' \donttest{
##'  library(competition)
##'  comp_mat = prep(data = euca, gen = 'clone', repl = 'block', area = 'area', 
##'                       ind = 'tree', age = 'age', row = 'row', col = 'col', 
##'                       dist.col = 3, dist.row = 2, trait = 'mai', method = 'SK',
##'                       n.dec = 3, verbose = TRUE)
##' }


prep<- function(data, gen, repl, row, col, ind, trait, dist.row, dist.col, 
                     method, area = NULL, age = NULL, n.dec = 2, 
                     verbose = FALSE){
  
  # Messages and warnings
  stopifnot("The number of individuals must be the product of no. rows * no. columns" = length(unique(data[,ind])) == 
              length(unique(data[,row])) * length(unique(data[,col])))
  stopifnot("Please, choose between the available methods ('MU', 'CC' or 'SK')" = method %in% c('MU', 'CC', 'SK'))
  
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
          n_dd = n_dd + 1
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
        gen = x[which(x$col == x$col[i] & x$row == x$row[[i]]),'trat'],
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
    
    class(Z) = "comprep"
    
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
            n_dd = n_dd + 1
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
          gen = q[which(q$col == q$col[i] & q$row == q$row[[i]]),'trat'],
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
    
    Z$data = do.call(rbind, lapply(Z, function(x) x$data))
    Z$control = control
    
    class(Z) = 'comprep'
    
    return(Z)
  }
}



#' Print an object of class `comprep`
#'
#' Print a `comprep` object in the R console. Alternatively, export the objects to the 
#' working directory
#'
#'
#' @param object An object of class `comprep`
#' @param category A string indicating which object to print. Options are "all" for
#' printing all objects, "data" for printing the data that will be used in the model, 
#' "matrix" for printing the competition matrix, "check" for printing the `neigh_check`
#' dataframe, and "CIF" for printing the competition intensity factor. 
#' @param age A string indicating if objects should be printed per age. Options 
#' are "all" for printing all ages, or the name of the specific age. Defaults to 
#' "all", which also serves when data has a single age.
#' 
#' @method print comprep
#' 
#' 
#' @importFrom data.table data.table
#' @importFrom utils write.csv
#' 
#' 
#' @export
#' 
#' @examples
#'\donttest{
#' library(competition)
#' comp_mat = prep(data = euca, gen = 'clone', repl = 'block', area = 'area', 
#'                       ind = 'tree', age = 'age', row = 'row', col = 'col', 
#'                       dist.col = 3, dist.row = 2, trait = 'mai', method = 'SK',
#'                       n.dec = 3, verbose = TRUE)
#' print(comp_mat, category = 'data', age = '6y')
#' print(comp_mat, category = 'matrix', age = 'all')
#' print(comp_mat, category = 'check', age = '3y')
#' print(comp_mat, category = 'CIF', age = 'all')
#' 
#' # Note that the ages are labeled as "3y" and "6y" in the example dataset 
#' }
#'


print.comprep = function(object, category = 'matrix', age = 'all'){
  
  stopifnot("The object must be of class 'comprep'" = class(object) == 'comprep')
  stopifnot("'category' should be of size 1" = length(category) == 1)
  stopifnot("'age' should be of size 1" = length(age) == 1)
  stopifnot("'age' does not exist" = age %in% c('all',gsub('Age_','',names(object)[grep("Age_*.", names(object))])))
  
  if(object$control[,6] == 0) age = 'all'
  
  # Data set
  
  if(category == 'all' | category == 'data'){
    
    if(age == 'all'){
      
      cat("\n","===> Data (competition matrix + user-provided data set)", "\n")
      print(data.table::data.table(object$data))
      
    }else{
      
      cat("\n","===> Data (competition matrix + user-provided data set), Age:", age, '\n')
      print(data.table::data.table(droplevels(object$data[which(object$data[,colnames(object$control)[6]] == age),])))
      
    }
  }
  
  # Competition matrix
  
  if(category == 'all' | category == 'matrix'){
    
    if(age == 'all'){
      
      cat("\n","===> Competition matrix", '\n')
      
      if(object$control[,6] == 0){
        
        print(object$Z)
        
      }else{
        
        for(i in gsub('Age_','',names(object)[grep("Age_*.", names(object))])){
          
          cat("\n","==> Age:", i, '\n')
          
          print(object[[grep(i, names(object))]][['Z']])
          
        }
      }
      
    }else{
      
      cat("\n","===> Competition matrix, Age:", age, '\n')
      print(object[[grep(age, names(object))]][['Z']])
      
    }
  }
  
  # Neighbourhood check
  
  if(category == 'all' | category == 'check'){
    
    if(age == 'all'){
      
      cat("\n","===> Neighbourhood check", "\n")
      
      if(object$control[,6] == 0){
        
      print(data.table::data.table(object$neigh_check))
        
      }else{
        
        for(i in gsub('Age_','',names(object)[grep("Age_*.", names(object))])){
          
          cat("\n","==> Age:", i, '\n')
          
          print(data.table(object[[grep(i, names(object))]][['neigh_check']]))
          
        }
      }
      
    }else{
      
      cat("\n","===> Neighbourhood check, Age:", age, '\n')
      print(data.table(object[[grep(age, names(object))]][['neigh_check']]))
      
    }
  }
  
  # CIF
  
  if(category == 'all' | category == 'CIF'){
    
    if(age == 'all'){
      
      cat("\n","===> Competition intensity factor", "\n")
      
      if(object$control[,6] == 0){
        
        print(object$CIF)
        
      }else{
        
        for(i in gsub('Age_','',names(object)[grep("Age_*.", names(object))])){
          
          cat("\n","==> Age:", i, '\n')
          
          print(object[[grep(i, names(object))]][['CIF']])
          
        }
      }
      
    }else{
      
      cat("\n","===> Competition intensity factor, Age:", age, '\n')
      print(object[[grep(age, names(object))]][['CIF']])
      
    }
  }

}


#' Print an object of class `comprep`
#'
#' Print a `comprep` object in the R console. 
#'
#'
#' @param object An object of class `comprep`
#' @param category A string indicating which object to print. Options are "all" for
#' printing all objects, "data" for printing the data that will be used in the model, 
#' "matrix" for printing the competition matrix, "check" for printing the `neigh_check`
#' dataframe, and "CIF" for printing the competition intensity factor. 
#' @param age A string indicating if objects should be printed per age. Options 
#' are "all" for printing all ages, or the name of the specific age. Defaults to 
#' "all", which also serves when data has a single age.
#' 
#' @method print comprep
#' 
#' 
#' @importFrom data.table data.table 
#' 
#' @export
#' 
#' @examples
#'\donttest{
#' library(competition)
#' comp_mat = prep(data = euca, gen = 'clone', repl = 'block', area = 'area', 
#'                       ind = 'tree', age = 'age', row = 'row', col = 'col', 
#'                       dist.col = 3, dist.row = 2, trait = 'mai', method = 'SK',
#'                       n.dec = 3, verbose = TRUE)
#' print(comp_mat, category = 'data', age = '6y')
#' print(comp_mat, category = 'matrix', age = 'all')
#' print(comp_mat, category = 'check', age = '3y')
#' print(comp_mat, category = 'CIF', age = 'all')
#' 
#' # Note that the ages are labeled as "3y" and "6y" in the example dataset 
#' }
#'


print.comprep = function(object, category = 'matrix', age = 'all'){
  
  stopifnot("The object must be of class 'comprep'" = class(object) == 'comprep')
  stopifnot("'category' should be of size 1" = length(category) == 1)
  stopifnot("'age' should be of size 1" = length(age) == 1)
  stopifnot("'age' does not exist" = age %in% c('all',gsub('Age_','',names(object)[grep("Age_*.", names(object))])))
  
  if(object$control[,6] == 0) age = 'all'
  
  # Data set
  
  if(category == 'all' | category == 'data'){
    
    if(age == 'all'){
      
      cat("\n","===> Data (competition matrix + user-provided data set)", "\n")
      print(data.table::data.table(object$data))
      
    }else{
      
      cat("\n","===> Data (competition matrix + user-provided data set), Age:", age, '\n')
      print(data.table::data.table(droplevels(object$data[which(object$data[,colnames(object$control)[6]] == age),])))
      
    }
  }
  
  # Competition matrix
  
  if(category == 'all' | category == 'matrix'){
    
    if(age == 'all'){
      
      cat("\n","===> Competition matrix", '\n')
      
      if(object$control[,6] == 0){
        
        print(object$Z)
        
      }else{
        
        for(i in gsub('Age_','',names(object)[grep("Age_*.", names(object))])){
          
          cat("\n","==> Age:", i, '\n')
          
          print(object[[grep(i, names(object))]][['Z']])
          
        }
      }
      
    }else{
      
      cat("\n","===> Competition matrix, Age:", age, '\n')
      print(object[[grep(age, names(object))]][['Z']])
      
    }
  }
  
  # Neighbourhood check
  
  if(category == 'all' | category == 'check'){
    
    if(age == 'all'){
      
      cat("\n","===> Neighbourhood check", "\n")
      
      if(object$control[,6] == 0){
        
        print(data.table::data.table(object$neigh_check))
        
      }else{
        
        for(i in gsub('Age_','',names(object)[grep("Age_*.", names(object))])){
          
          cat("\n","==> Age:", i, '\n')
          
          print(data.table(object[[grep(i, names(object))]][['neigh_check']]))
          
        }
      }
      
    }else{
      
      cat("\n","===> Neighbourhood check, Age:", age, '\n')
      print(data.table(object[[grep(age, names(object))]][['neigh_check']]))
      
    }
  }
  
  # CIF
  
  if(category == 'all' | category == 'CIF'){
    
    if(age == 'all'){
      
      cat("\n","===> Competition intensity factor", "\n")
      
      if(object$control[,6] == 0){
        
        print(object$CIF)
        
      }else{
        
        for(i in gsub('Age_','',names(object)[grep("Age_*.", names(object))])){
          
          cat("\n","==> Age:", i, '\n')
          
          print(object[[grep(i, names(object))]][['CIF']])
          
        }
      }
      
    }else{
      
      cat("\n","===> Competition intensity factor, Age:", age, '\n')
      print(object[[grep(age, names(object))]][['CIF']])
      
    }
  }
  
}


#' Plots for data overview
#'
#' The function generates two types of plots: i) a heatmap representing the grid. 
#' The cells are filled according to the phenotype value of each plot; and ii) 
#' boxplots depicting the performance of each selection candidate
#'
#'
#' @param object An object of class `comprep`
#' @param category A string indicating which plot to build. Options are "heatmap" for
#' plotting the field grid, or "boxplot" for plotting the boxplots of phenotypic 
#' performance.
#' 
#' @method plot comprep
#' 
#' @details
#' All plots are built using the [ggplot2::ggplot()] library, so they are all 
#' customizable using "+ ggfun"
#' 
#' 
#' @importFrom ggplot2 ggplot
#' 
#' @export
#' 
#' @examples
#'\donttest{
#' library(competition)
#' comp_mat = prep(data = euca, gen = 'clone', repl = 'block', area = 'area', 
#'                       ind = 'tree', age = 'age', row = 'row', col = 'col', 
#'                       dist.col = 3, dist.row = 2, trait = 'mai', method = 'SK',
#'                       n.dec = 3, verbose = TRUE)
#'                       
#' 
#' }
#'

plot.comprep = function(object, category){
 
  stopifnot("The object must be of class 'comprep'" = class(object) == 'comprep')
  stopifnot("'category' should be of size 1" = length(category) == 1)
  stopifnot("Please, choose between the available categories ('heatmap' or 'boxplot')" = category %in% c('heatmap', 'boxplot'))
  
  
  if(category == 'heatmap'){
    
    if(object$control[,6] > 0){
      if(object$control[,7] > 0){
        
        dat = as.data.frame(object$data)[,names(object$control)[c(2,4,5,6,7,1)]]
        
        facet.label.col = paste(names(object$control)[7], unique(dat[,names(object$control)[7]]))
        names(facet.label.col) = unique(dat[,names(object$control)[7]])
        facet.label.row =  paste(names(object$control)[6], unique(dat[,names(object$control)[6]]))
        names(facet.label.row) = unique(dat[,names(object$control)[6]])
        
        colnames(dat) = c('gen', 'row', 'col', 'age', 'area', 'y')
        
        ggplot(data = dat, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$y)) + 
          geom_tile(color = 'black') + 
          facet_grid(rows = vars(.data$age), cols = vars(.data$area), 
                     scales = 'free_x', 
                     labeller = labeller(.cols = facet.label.col,
                                         .rows = facet.label.row)) +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = 'Y') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))
        
        
      } else {
        
        dat = as.data.frame(object$data)[,names(object$control)[c(2,4,5,6,1)]]
        
        facet.label.row =  paste(names(object$control)[6], unique(dat[,names(object$control)[6]]))
        names(facet.label.row) = unique(dat[,names(object$control)[6]])
        
        colnames(dat) = c('gen', 'row', 'col', 'age', 'y')
        
        ggplot(data = dat, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$y)) + 
          geom_tile(color = 'black') + 
          facet_grid(rows = vars(.data$age),
                     scales = 'free_x', 
                     labeller = labeller(.rows = facet.label.row)) +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = 'Y') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))
        
      }
    } else {
      if(object$control[,7] > 0){
        
        dat = as.data.frame(object$data)[,names(object$control)[c(2,4,5,7,1)]]
        
        facet.label.col = paste(names(object$control)[7], unique(dat[,names(object$control)[7]]))
        names(facet.label.col) = unique(dat[,names(object$control)[7]])
        
        colnames(dat) = c('gen', 'row', 'col', 'area', 'y')
        
        ggplot(data = dat, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$y)) + 
          geom_tile(color = 'black') + 
          facet_grid(cols = vars(.data$area), 
                     scales = 'free_x', 
                     labeller = labeller(.cols = facet.label.col)) +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = 'Y') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))

        
      } else {
        
        dat = as.data.frame(object$data)[,names(object$control)[c(2,4,5,1)]]

        colnames(dat) = c('gen', 'row', 'col', 'y')
        
        ggplot(data = dat, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$y)) + 
          geom_tile(color = 'black') + 
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = 'Y') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))

      }
    }
  }else if(category == 'boxplot'){

    if(object$control[,6] > 0){
      if(object$control[,7] > 0){
        
        dat = as.data.frame(object$data)[,names(object$control)[c(2,4,5,6,7,1)]]
        
        facet.label.col = paste(names(object$control)[7], unique(dat[,names(object$control)[7]]))
        names(facet.label.col) = unique(dat[,names(object$control)[7]])
        facet.label.row =  paste(names(object$control)[6], unique(dat[,names(object$control)[6]]))
        names(facet.label.row) = unique(dat[,names(object$control)[6]])
        
        colnames(dat) = c('gen', 'row', 'col', 'age', 'area', 'y')
        
        suppressWarnings({
          ggplot(data = dat, aes(x = .data$gen, y = .data$y)) + 
            geom_boxplot() + 
            theme_bw() + 
            theme(axis.text.x = element_text(angle = 90)) + 
            facet_grid(rows = vars(.data$age), cols = vars(.data$area), 
                       labeller = labeller(.cols = facet.label.col,
                                           .rows = facet.label.row)) + 
            labs(x = 'Genotype', y = 'Y')
        })

        
      } else {
        
        dat = as.data.frame(object$data)[,names(object$control)[c(2,4,5,6,1)]]
        
        facet.label.row =  paste(names(object$control)[6], unique(dat[,names(object$control)[6]]))
        names(facet.label.row) = unique(dat[,names(object$control)[6]])
        
        colnames(dat) = c('gen', 'row', 'col', 'age', 'y')
        
        suppressWarnings({
          ggplot(data = dat, aes(x = .data$gen, y = .data$y)) + 
            geom_boxplot() + 
            theme_bw() + 
            theme(axis.text.x = element_text(angle = 90)) + 
            facet_grid(rows = vars(.data$age), 
                       labeller = labeller(.rows = facet.label.row)) + 
            labs(x = 'Genotype', y = 'Y')
        })
        
      }
    } else {
      if(object$control[,7] > 0){
        
        dat = as.data.frame(object$data)[,names(object$control)[c(2,4,5,7,1)]]
        
        facet.label.col = paste(names(object$control)[7], unique(dat[,names(object$control)[7]]))
        names(facet.label.col) = unique(dat[,names(object$control)[7]])
        
        colnames(dat) = c('gen', 'row', 'col', 'area', 'y')
        
        suppressWarnings({
          ggplot(data = dat, aes(x = .data$gen, y = .data$y)) + 
            geom_boxplot() + 
            theme_bw() + 
            theme(axis.text.x = element_text(angle = 90)) + 
            facet_grid(cols = vars(.data$area), 
                       labeller = labeller(.cols = facet.label.col)) + 
            labs(x = 'Genotype', y = 'Y')
        })
        
        
      } else {
        
        dat = as.data.frame(object$data)[,names(object$control)[c(2,4,5,1)]]
        
        colnames(dat) = c('gen', 'row', 'col', 'y')
        
        suppressWarnings({
          ggplot(data = dat, aes(x = .data$gen, y = .data$y)) + 
            geom_boxplot() + 
            theme_bw() + 
            theme(axis.text.x = element_text(angle = 90)) + 
            labs(x = 'Genotype', y = 'Y')
        })
        
      }
    }

  }

}


#Field
# x$trat_cod = ifelse(is.na(x$trait), paste0(x$trat,'*'), paste(x$trat))
# field = stats::reshape(data = x[,c('row', 'col', 'trat_cod')], direction = 'wide', 
#                        v.names = 'trat_cod', timevar = 'col', idvar = 'row')
# rownames(field) = field$row; field = field[,-1]
# colnames(field) = sub('trat_cod.', '', colnames(field))