input = prep.out$data
control = attr(prep.out, 'control')

input = input[order(input[,names(control)[7]], input[,names(control)[6]],
                    input[,names(control)[4]], input[,names(control)[5]]),]

scm1 = asreml::asreml(fixed = fixed,
                     random = stats::as.formula(
                       paste(paste(as.character(random), collapse = ""), 
                             paste0('+ str(~', colnames(control)[2],
                                    '+ grp(g1), ~corh(2):id(',
                                    colnames(control)[2],')) + str(~',
                                    colnames(control)[2],':',
                                    colnames(control)[6],'+grp(g1):',
                                    colnames(control)[6], ', ~diag(2):id(',
                                    colnames(control)[6], '):id(',
                                    colnames(control)[2], '))'))
                     ),
                     group = list(g1 = 1:control[,2]),
                     residual = stats::as.formula(paste0(
                       '~dsum(~id(', colnames(control)[6],'):ar1(',
                       colnames(control)[4],'):ar1(',
                       colnames(control)[5],')|',
                       colnames(control)[7],')')
                     ),
                     data = input)


scm2 = asreml::asreml(fixed = fixed,
                      random = stats::as.formula(
                        paste(paste(as.character(random), collapse = ""), 
                              paste0('+ str(~', colnames(control)[2],
                                     '+ grp(g1), ~corh(2):id(',
                                     colnames(control)[2],')) + str(~',
                                     colnames(control)[2],':',
                                     colnames(control)[6],'+grp(g1):',
                                     colnames(control)[6], ', ~diag(2):id(',
                                     colnames(control)[6], '):id(',
                                     colnames(control)[2], '))'))
                      ),
                      group = list(g1 = 1:control[,2]),
                      residual = stats::as.formula(paste0(
                        '~dsum(~diag(', colnames(control)[6],'):ar1(',
                        colnames(control)[4],'):ar1(',
                        colnames(control)[5],')|',
                        colnames(control)[7],')')
                      ),
                      data = input)


scm3 = asreml::asreml(fixed = fixed,
                      random = stats::as.formula(
                        paste(paste(as.character(random), collapse = ""), 
                              paste0('+ str(~', colnames(control)[2],
                                     '+ grp(g1), ~corh(2):id(',
                                     colnames(control)[2],')) + str(~',
                                     colnames(control)[2],':',
                                     colnames(control)[6],'+grp(g1):',
                                     colnames(control)[6], ', ~diag(2):id(',
                                     colnames(control)[6], '):id(',
                                     colnames(control)[2], '))'))
                      ),
                      group = list(g1 = 1:control[,2]),
                      residual = stats::as.formula(paste0(
                        '~dsum(~ar1(', colnames(control)[6],'):ar1(',
                        colnames(control)[4],'):ar1(',
                        colnames(control)[5],')|',
                        colnames(control)[7],')')
                      ),
                      data = input)

scm4 = asreml::asreml(fixed = fixed,
                      random = stats::as.formula(
                        paste(paste(as.character(random), collapse = ""), 
                              paste0('+ str(~', colnames(control)[2],
                                     '+ grp(g1), ~corh(2):id(',
                                     colnames(control)[2],')) + str(~',
                                     colnames(control)[2],':',
                                     colnames(control)[6],'+grp(g1):',
                                     colnames(control)[6], ', ~diag(2):id(',
                                     colnames(control)[6], '):id(',
                                     colnames(control)[2], '))'))
                      ),
                      group = list(g1 = 1:control[,2]),
                      residual = stats::as.formula(paste0(
                        '~dsum(~ar1h(', colnames(control)[6],'):ar1(',
                        colnames(control)[4],'):ar1(',
                        colnames(control)[5],')|',
                        colnames(control)[7],')')
                      ),
                      data = input, maxit = 50)


data.frame(
  model = c("ID", "DIAG", "AR1", "AR1H"),
  AIC = c(summary(scm1)$aic, summary(scm2)$aic, summary(scm3)$aic, summary(scm4)$aic),
  BIC = c(summary(scm1)$bic, summary(scm2)$bic, summary(scm3)$bic, summary(scm4)$bic),
  LogL = c(summary(scm1)$logl, summary(scm2)$logl, summary(scm3)$logl, summary(scm4)$logl)
)

summary(scm1)$varcomp
summary(scm2)$varcomp
summary(scm3)$varcomp
summary(scm4)$varcomp




