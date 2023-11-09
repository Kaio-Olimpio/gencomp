##'
##'
##'
##'
##'
##'
##'
##'
##' @import asreml
##' 
##' @export


gscmod = function(gscprep.out, fixed, random = ~1, cor = TRUE, maxit = 50, ...) {
  
  requireNamespace('asreml')
  
  gscprep.out <<- gscprep.out
  
  input = gscprep.out$data
  
  if(gscprep.out$control[,7] > 0) ### Multi-areas ------------
    {
      if(gscprep.out$control[,6] > 0) ### Multi-ages --------
      {
        if(cor){
          input = input[order(input$area, input$age, input$row, input$col),]
          
          scm = asreml::asreml(fixed = fixed,
                       random = as.formula(
                         paste(paste(as.character(random), collapse = ""), 
                               paste0('+ str(~', colnames(gscprep.out$control)[2],
                                      '+ grp(g1), ~corh(2):id(',
                                      colnames(gscprep.out$control)[2],')) + str(~',
                                      colnames(gscprep.out$control)[2],':',
                                      colnames(gscprep.out$control)[6],'+grp(g1):',
                                      colnames(gscprep.out$control)[6], ', ~corh(2):id(',
                                      colnames(gscprep.out$control)[6], '):id(',
                                      colnames(gscprep.out$control)[2], '))'))
                       ),
                       group = list(g1 = 1:gscprep.out$control[,2]),
                       residual = as.formula(paste0(
                         '~dsum(~ar1h(', colnames(gscprep.out$control)[6],'):ar1(',
                         colnames(gscprep.out$control)[4],'):ar1(',
                         colnames(gscprep.out$control)[5],')|',
                         colnames(gscprep.out$control)[7],')')
                       ),
                       maxit = maxit,
                       data = input, ...)

        }else{
          input = input[order(input$area, input$age, input$row, input$col),]
          
          scm = asreml::asreml(fixed = fixed,
                       random = as.formula(
                         paste(paste(as.character(random), collapse = ""), 
                               paste0('+', colnames(gscprep.out$control)[2],
                                      '+ grp(g1) +',
                                      colnames(gscprep.out$control)[2],':',
                                      colnames(gscprep.out$control)[6],'+grp(g1):',
                                      colnames(gscprep.out$control)[6]))
                       ),
                       group = list(g1 = 1:gscprep.out$control[,2]),
                       residual = as.formula(paste0(
                         '~dsum(~ar1h(', colnames(gscprep.out$control)[6],'):ar1(',
                         colnames(gscprep.out$control)[4],'):ar1(',
                         colnames(gscprep.out$control)[5],')|',
                         colnames(gscprep.out$control)[7],')')
                       ),
                       maxit = maxit,
                       data = input, ...)

        }
      }else ### Single age ---------------
      {
        if(cor){
          scm = asreml::asreml(fixed = fixed,
                       random = as.formula(
                         paste(paste(as.character(random), collapse = ""), 
                               paste0('+ str(~', colnames(gscprep.out$control)[2],
                                      '+ grp(g1), ~corh(2):id(',
                                      colnames(gscprep.out$control)[2],'))'))
                       ),
                       group = list(g1 = 1:gscprep.out$control[,2]),
                       residual = as.formula(paste0(
                         '~dsum(~ar1(', colnames(gscprep.out$control)[4],'):ar1(',
                         colnames(gscprep.out$control)[5],')|',
                         colnames(gscprep.out$control)[7],')')
                       ),
                       maxit = maxit,
                       data = input, ...)
        }else{
          scm = asreml::asreml(fixed = fixed,
                       random = as.formula(
                         paste(paste(as.character(random), collapse = ""), 
                               paste0('+', colnames(gscprep.out$control)[2],
                                      '+ grp(g1)'))
                       ),
                       group = list(g1 = 1:gscprep.out$control[,2]),
                       residual = as.formula(paste0(
                         '~dsum(~ar1(', colnames(gscprep.out$control)[4],'):ar1(',
                         colnames(gscprep.out$control)[5],')|',
                         colnames(gscprep.out$control)[7],')')
                       ),
                       maxit = maxit,
                       data = input, ...)
        }
      } 
    }else ### Single area -----------
      {
        if(gscprep.out$control[,6] > 0) ### Multi-ages --------
        {
          if(cor){
            input = input[order(input$age, input$row, input$col),]
            scm = asreml::asreml(fixed = fixed,
                         random = as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+ str(~', colnames(gscprep.out$control)[2],
                                        '+ grp(g1), ~corh(2):id(',
                                        colnames(gscprep.out$control)[2],')) + str(~',
                                        colnames(gscprep.out$control)[2],':',
                                        colnames(gscprep.out$control)[6],'+grp(g1):',
                                        colnames(gscprep.out$control)[6], ', ~corh(2):id(',
                                        colnames(gscprep.out$control)[6], '):id(',
                                        colnames(gscprep.out$control)[2], '))'))
                         ),
                         group = list(g1 = 1:gscprep.out$control[,2]),
                         residual = as.formula(paste0(
                           '~ar1h(', colnames(gscprep.out$control)[6],'):ar1(',
                           colnames(gscprep.out$control)[4],'):ar1(',
                           colnames(gscprep.out$control)[5],')')
                         ),
                         maxit = maxit,
                         data = input, ...)
            
          }else{
            input = input[order(input$age, input$row, input$col),]
            scm = asreml::asreml(fixed = fixed,
                         random = as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(gscprep.out$control)[2],
                                        '+ grp(g1) +',
                                        colnames(gscprep.out$control)[2],':',
                                        colnames(gscprep.out$control)[6],'+grp(g1):',
                                        colnames(gscprep.out$control)[6]))
                         ),
                         group = list(g1 = 1:gscprep.out$control[,2]),
                         residual = as.formula(paste0(
                           '~ar1h(', colnames(gscprep.out$control)[6],'):ar1(',
                           colnames(gscprep.out$control)[4],'):ar1(',
                           colnames(gscprep.out$control)[5],')')
                         ),
                         maxit = maxit,
                         data = input, ...)
          }
        }else ### Single age ---------------
        {
          if(cor){
            scm = asreml::asreml(fixed = fixed,
                         random = as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+ str(~', colnames(gscprep.out$control)[2],
                                        '+ grp(g1), ~corh(2):id(',
                                        colnames(gscprep.out$control)[2],'))'))
                         ),
                         group = list(g1 = 1:gscprep.out$control[,2]),
                         residual = as.formula(paste0(
                           '~ar1(', colnames(gscprep.out$control)[4],'):ar1(',
                           colnames(gscprep.out$control)[5],')')
                         ),
                         maxit = maxit,
                         data = input, ...)
          }else{
            scm = asreml::asreml(fixed = fixed,
                         random = as.formula(
                           paste(paste(as.character(random), collapse = ""), 
                                 paste0('+', colnames(gscprep.out$control)[2],
                                        '+ grp(g1)'))
                         ),
                         group = list(g1 = 1:gscprep.out$control[,2]),
                         residual = as.formula(paste0(
                           '~ar1(', colnames(gscprep.out$control)[4],'):ar1(',
                           colnames(gscprep.out$control)[5],')')
                         ),
                         maxit = maxit,
                         data = input, ...)
          }
        }
  }
  
  if(!scm$converge) warning("The model did not converge")
  
  remove(gscprep.out, envir = .GlobalEnv)
  
  return(scm)
}


