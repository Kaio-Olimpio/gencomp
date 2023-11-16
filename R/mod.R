##'
##'
##'
##'
##'
##'
##'
##'
##' @import asreml
##' @importFrom stats as.formula
##' 
##' @export
##' 

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


