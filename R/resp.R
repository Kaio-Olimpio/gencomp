##'
##'
##'
##'
##'
##'
##' @import ggplot2
##' @importFrom rlang .data
##' @importFrom stats reorder
##' @importFrom ggpubr ggarrange
##' 
##' @export
##' 
##' 


comp.resp = function(prep.out, model, weight.tgv = FALSE) {
  
  requireNamespace('ggplot2')
  
  prep.out <<- prep.out
  model = model
  output = list()
  
  # Variance components ---------------------
  varcomp = summary(model)$varcomp
  
  ## Dealing with the names --------------------
  if(prep.out$control[,6] > 0){
    varcomp = varcomp[-grep('!R', rownames(varcomp)),]
    rownames(varcomp)[
      rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
                                                    grepl('_1', rownames(varcomp)) &
                                                    grepl(names(prep.out$control)[6],
                                                          rownames(varcomp))),])
    ] = paste("DGE", names(prep.out$control)[6], sep = ':')
    rownames(varcomp)[
      rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
                                                    grepl('_2', rownames(varcomp)) &
                                                    grepl(names(prep.out$control)[6],
                                                          rownames(varcomp))),])
    ] = paste("IGE", names(prep.out$control)[6], sep = ':')
    # rownames(varcomp)[
    #   rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
    #                                                 grepl('cor', rownames(varcomp)) &
    #                                                 grepl(names(prep.out$control)[6],
    #                                                       rownames(varcomp))),])
    # ] = paste("cor(IGE_DGE)", names(prep.out$control)[6], sep = ':')
    
    if(prep.out$control[,7] > 0){
      rownames(varcomp)[
        rownames(varcomp) %in% rownames(varcomp[which(grepl(paste(names(prep.out$control)[6], 'cor', sep='!'), 
                                                            rownames(varcomp))),])
      ] = paste(paste0("R=autocor(", 
                       paste(names(prep.out$control[7]),
                             1:as.numeric(prep.out$control[7]), sep = '_'), ')'),
                names(prep.out$control)[6], sep = '!')
      
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl(paste(names(prep.out$control)[6], 
                                             '_', sep=''),
                                       rownames(varcomp))),])
      ] = paste0("R=", rownames(varcomp[which(grepl(paste(names(prep.out$control)[6], 
                                                          '_', sep=''),
                                                    rownames(varcomp))),]))
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('row!',
                                       rownames(varcomp))),])
      ] = paste0('R=autocor(row):', 
                 sub('!.*', '', 
                     rownames(varcomp[which(grepl('row',rownames(varcomp))),])))
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('col',
                                       rownames(varcomp))),])
      ] = paste0('R=autocor(col):', 
                 sub('!.*', '', 
                     rownames(varcomp[which(grepl('col',rownames(varcomp))),])))
      
    }else{
      rownames(varcomp)[
        rownames(varcomp) == rownames(varcomp[which(grepl(paste(names(prep.out$control)[6], 'cor', sep='!'), 
                                                          rownames(varcomp))),])
      ] = paste0("R=autocor(", names(prep.out$control)[6],')')
      
      rownames(varcomp)[
        rownames(varcomp) == rownames(varcomp[which(grepl(paste(names(prep.out$control)[6], '_', sep=''), 
                                                          rownames(varcomp))),] )
      ] = paste0("R=", paste(names(prep.out$control)[6],
                             1:as.numeric(prep.out$control[6]), sep = '_'))
      
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('row!',
                                       rownames(varcomp))),])
      ] = 'R=autocor(row)'
      
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('col',
                                       rownames(varcomp))),])
      ] = 'R=autocor(col)'
      
    }
  }else{
    
    if(prep.out$control[,7] > 0){
      rownames(varcomp)[
        rownames(varcomp) %in% rownames(varcomp[which(grepl('!R', 
                                                            rownames(varcomp))),])
      ] = paste0("R=", 
                 sub('!.*', '',
                     rownames(varcomp[which(grepl('col',rownames(varcomp))),])))
      
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('row!',
                                       rownames(varcomp))),])
      ] = paste0('R=autocor(row):', 
                 sub('!.*', '', 
                     rownames(varcomp[which(grepl('row',rownames(varcomp))),])))
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('col',
                                       rownames(varcomp))),])
      ] = paste0('R=autocor(col):', 
                 sub('!.*', '', 
                     rownames(varcomp[which(grepl('col',rownames(varcomp))),])))
      
    }else{
      rownames(varcomp)[
        rownames(varcomp) %in% rownames(varcomp[which(grepl('!R', 
                                                            rownames(varcomp))),])
      ] = 'R'
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('row!',
                                       rownames(varcomp))),])
      ] = 'R=autocor(row):'
      rownames(varcomp)[
        rownames(varcomp) %in% 
          rownames(varcomp[which(grepl('col',
                                       rownames(varcomp))),])
      ] = 'R=autocor(col):'
    }
    
  }
  
  rownames(varcomp)[
    rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
                                                  grepl('_1', rownames(varcomp))),])
  ] = 'DGE'
  rownames(varcomp)[
    rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
                                                  grepl('_2', rownames(varcomp))),])
  ] = 'IGE'
  rownames(varcomp)[
    rownames(varcomp) == rownames(varcomp[which(grepl('grp', rownames(varcomp)) &
                                                  grepl('cor', rownames(varcomp))),])
  ] = "cor(IGE_DGE)"
  
  output$varcomp = varcomp
  
  
  # BLUPs --------------------
  blup = as.data.frame(summary(model, coef = TRUE)$coef.random)
  
  if(prep.out$control[,6] > 0){
    ## Main effects --------------
    DGE = blup[which(grepl(names(prep.out$control)[2], rownames(blup)) & 
                       !grepl(names(prep.out$control)[6], rownames(blup))), -3]
    DGE[,names(prep.out$control)[2]] = gsub(paste0(names(prep.out$control)[2],'_'), 
                                               '', rownames(DGE), fixed = T)
    rownames(DGE) = NULL
    DGE$rel.DGE = 1-(DGE$std.error^2/varcomp['DGE','component'])
    DGE = DGE[,c(3,1,2,4)]; colnames(DGE)[c(2,3)] = c('DGE', 'se.DGE')
    
    IGE = blup[which(grepl('grp', rownames(blup)) & 
                       !grepl(names(prep.out$control)[6], rownames(blup))), -3]
    IGE[,names(prep.out$control)[2]] = regmatches(rownames(IGE), 
                                                  m = regexpr(paste(DGE$clone, collapse = "|"),
                                                              rownames(IGE)))
    rownames(IGE) = NULL
    IGE$rel.IGE = 1-(IGE$std.error^2/varcomp['IGE','component'])
    IGE = IGE[,c(3,1,2,4)]; colnames(IGE)[c(2,3)] = c('IGE', 'se.IGE')
    
    IGE$class = ifelse(IGE$IGE > mean(IGE$IGE) + sd(IGE$IGE), 'Sensitive', 
                       ifelse(IGE$IGE < mean(IGE$IGE) - sd(IGE$IGE), 'Aggressive', 
                              'Homeostatic'))

    main = merge(DGE, IGE, by = names(prep.out$control)[2])
    
    if(isFALSE(weight.tgv)){
      main$TGV = main$DGE + mean(do.call(c,lapply(Filter(function(x) is.list(x), 
                                              prep.out), function(x){x$CIF}))) *
        main$IGE
    }else if(is.numeric(weight.tgv)){
      main$TGV = weight.tgv[1] * main$DGE + 
        mean(do.call(c,lapply(Filter(function(x) is.list(x),
                                     prep.out), function(x){x$CIF}))) *
        main$IGE * weight.tgv[2]
    }else if(isTRUE(weight.tgv)){
      main$TGV = main$rel.DGE * main$DGE + 
        mean(do.call(c,lapply(Filter(function(x) is.list(x),
                                     prep.out), function(x){x$CIF}))) *
        main$IGE * main$rel.IGE
    }
    
    ## Interaction effects ----------------
    DGE.int = blup[which(grepl(names(prep.out$control)[2], rownames(blup)) & 
                           grepl(names(prep.out$control)[6], rownames(blup))), -3]
    DGE.int[,names(prep.out$control)[2]] = gsub(paste0(names(prep.out$control)[2],'_'),'', 
                                                gsub(':.*', '', rownames(DGE.int)), fixed = T)
    DGE.int[,names(prep.out$control)[6]] = gsub(paste0(names(prep.out$control)[6], '_'),'', 
                                                gsub('.*:', '', rownames(DGE.int)), fixed = T)
    rownames(DGE.int) = NULL
    DGE.int$rel.DGE = 1-(DGE.int$std.error^2/varcomp['DGE','component'])
    DGE.int = DGE.int[,c(3,4,1,2,5)]; colnames(DGE.int)[c(3,4)] = c('DGE', 'se.DGE')
    DGE.int = merge(DGE.int, DGE[,c(names(prep.out$control)[2], 'DGE')],
                    by = names(prep.out$control)[2])
    DGE.int$DGE = DGE.int$DGE.x + DGE.int$DGE.y
    DGE.int = DGE.int[,c(1, 2, 7, 4, 5)]
    
    IGE.int = blup[which(grepl('grp', rownames(blup)) & 
                           grepl(names(prep.out$control)[6], rownames(blup))), -3]
    IGE.int[,names(prep.out$control)[2]] = regmatches(rownames(IGE.int), 
                                                      m = regexpr(paste(DGE$clone, collapse = "|"),
                                                                  rownames(IGE.int)))
    IGE.int[,names(prep.out$control)[6]] = gsub(paste0(names(prep.out$control)[6], '_'),'', 
                                                   gsub('.*:', '', rownames(IGE.int)), fixed = T)
    rownames(IGE.int) = NULL
    IGE.int$rel.IGE = 1-(IGE.int$std.error^2/varcomp['IGE','component'])
    IGE.int = IGE.int[,c(3,4,1,2,5)]; colnames(IGE.int)[c(3,4)] = c('IGE', 'se.IGE')
    IGE.int = merge(IGE.int, IGE[,c(names(prep.out$control)[2], 'IGE')],
                    by = names(prep.out$control)[2])
    IGE.int$IGE = IGE.int$IGE.x + IGE.int$IGE.y
    IGE.int = IGE.int[,c(1, 2, 7, 4, 5)]
    IGE.int$class = do.call(c, lapply(split(IGE.int, IGE.int[,2]), function(x){
      ifelse(x$IGE > mean(x$IGE) + sd(x$IGE), 'Sensitive', 
             ifelse(x$IGE < mean(x$IGE) - sd(x$IGE), 'Aggressive', 
                    'Homeostatic'))
    }))
    
    within = merge(DGE.int, IGE.int, 
                   by = c(names(prep.out$control)[2],
                          names(prep.out$control)[6]))

    if(isFALSE(weight.tgv)){
      within =  merge(within, do.call(rbind, lapply(split(within, within[,2]), function(x){
        data.frame(
          gen = x[,1], 
          age = x[,2],
          TGV = x$DGE + prep.out[[grep(unique(x[,2]),names(prep.out))]]$CIF *
            x$IGE
        )
      })), by.x = c(names(prep.out$control)[2], names(prep.out$control)[6]), 
      by.y = c('gen','age'))
    }else if(is.numeric(weight.tgv)){
      within =  merge(within, do.call(rbind, lapply(split(within, within[,2]), function(x){
        data.frame(
          gen = x[,1], 
          age = x[,2],
          TGV =  weight.tgv[1] * x$DGE + 
            prep.out[[grep(unique(x[,2]),names(prep.out))]]$CIF *
            x$IGE* weight.tgv[2]
        )
      })), by.x = c(names(prep.out$control)[2], names(prep.out$control)[6]), 
      by.y = c('gen','age'))
    }else if(isTRUE(weight.tgv)){
      within =  merge(within, do.call(rbind, lapply(split(within, within[,2]), function(x){
        data.frame(
          gen = x[,1], 
          age = x[,2],
          TGV =  x$rel.DGE * x$DGE + 
            prep.out[[grep(unique(x[,2]),names(prep.out))]]$CIF *
            x$IGE* x$rel.IGE
        )
      })), by.x = c(names(prep.out$control)[2], names(prep.out$control)[6]), 
      by.y = c('gen','age'))
    }
    
    output$blups = list(main = main, within = within)
    
    ## Other random effects ----------------
    if(any(!grepl("grp\\(g1\\)", attr(model$formulae$random, 'term.labels')))) {
      coef.random = blup[which(!grepl('grp', rownames(blup)) & 
                                 !grepl(names(prep.out$control)[2], rownames(blup))), -3]
      output$blups$coef.random = coef.random
    }
  }else{
    
    ## Main effects --------------
    DGE = blup[which(grepl(names(prep.out$control)[2], rownames(blup)) & 
                       !grepl(names(prep.out$control)[6], rownames(blup))), -3]
    DGE[,names(prep.out$control)[2]] = gsub(paste0(names(prep.out$control)[2],'_'), 
                                               '', rownames(DGE), fixed = T)
    rownames(DGE) = NULL
    DGE$rel.DGE = 1-(DGE$std.error^2/varcomp['DGE','component'])
    DGE = DGE[,c(3,1,2,4)]; colnames(DGE)[c(2,3)] = c('DGE', 'se.DGE')
    
    IGE = blup[which(grepl('grp', rownames(blup)) & 
                       !grepl(names(prep.out$control)[6], rownames(blup))), -3]
    IGE[,names(prep.out$control)[2]] = regmatches(rownames(IGE), 
                                                  m = regexpr(paste(DGE$clone, collapse = "|"), rownames(IGE)))
    rownames(IGE) = NULL
    IGE$rel.IGE = 1-(IGE$std.error^2/varcomp['IGE','component'])
    IGE = IGE[,c(3,1,2,4)]; colnames(IGE)[c(2,3)] = c('IGE', 'se.IGE')
    
    IGE$class = ifelse(IGE$IGE > mean(IGE$IGE) + sd(IGE$IGE), 'Sensitive', 
                       ifelse(IGE$IGE < mean(IGE$IGE) - sd(IGE$IGE), 'Aggressive', 
                              'Homeostatic'))
    
    main = merge(DGE, IGE, by = names(prep.out$control)[2])    
    
    if(isFALSE(weight.tgv)){
      main$TGV = main$DGE + prep.out$CIF *  main$IGE
    }else if(is.numeric(weight.tgv)){
      main$TGV = weight.tgv[1] * main$DGE + prep.out$CIF *
        main$IGE * weight.tgv[2]
    }else if(isTRUE(weight.tgv)){
      main$TGV = main$rel.DGE * main$DGE + prep.out$CIF *
        main$IGE * main$rel.IGE
    }
    
    ## Other random effects ----------------
    if(any(!grepl("grp\\(g1\\)",attr(model$formulae$random, 'term.labels')))) {
      coef.random = blup[which(!grepl('grp', rownames(blup)) & 
                                 !grepl(names(prep.out$control)[2], rownames(blup))), -3]
      output$blups = list(main = main, coef.random = coef.random)
    }else{
      output$blups = main
    }
  }
  
  # Plots -----------------------------
  ## Main effects
  output$plots$main = list()
  
  ### Density
  den.ige.main = ggplot(data = main) + 
    geom_density(aes(x = .data$IGE, after_stat(density), fill = 'Homeostatic'), 
                 alpha = .8)
  d = ggplot2::ggplot_build(den.ige.main)$data[[1]]
  den.ige.main = den.ige.main + geom_area(data = subset(d, x > mean(IGE$IGE) + sd(IGE$IGE)), 
                aes(x = .data$x, y = .data$y, fill = 'Sensitive'), alpha = 0.8) + 
    geom_area(data = subset(d, x < mean(IGE$IGE) - sd(IGE$IGE)), 
              aes(x = .data$x, y = .data$y, fill = 'Aggressive'), alpha = .8) +
    geom_density(data = main, aes(x = .data$IGE, after_stat(density)), 
                 color = 'black', linewidth = 1.2, show.legend = F) +
    scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
    labs(x = "IGE", y = 'Density', fill = 'Class') + theme_minimal() +
    theme(legend.position = 'top')
  
  output$plots$main$IGE.density = den.ige.main
  
  ### DGE vs IGE
  plott = ggplot(data = main, aes(x = .data$IGE, y = .data$DGE)) + 
    geom_point(aes(fill = .data$class), colour = 'black', pch = 21, size = 2.2) +
    scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
    labs(fill = 'Class') + theme_minimal() +
    theme(legend.position = 'top') 
  
  output$plots$main$DGEvsIGE = plott
  
  ### DGE and IGE
  temp = data.frame(
    gen = rep(main[, 1], times = 2), 
    blup = c(main$DGE, main$IGE),
    se = c(main$se.DGE, main$se.IGE),
    rel = c(main$rel.DGE, main$rel.IGE),
    comp = rep(c('DGE', 'IGE'), each = prep.out$control[,2])
  )
  
  temp$gen = factor(temp$gen, levels = temp[order(temp$blup, decreasing = TRUE)[
    order(temp$blup, decreasing = TRUE) %in% which(temp$comp == 'DGE')], 1
    ])
  
  plott = ggplot(data = temp, aes(x = .data$gen, y = .data$blup)) + 
    facet_wrap(.~.data$comp, ncol = 1, scales = 'free_y') + 
    geom_segment(aes(x = .data$gen,
                     xend = .data$gen, y = 0, yend = .data$blup),
                 linewidth = 1) + 
    geom_point(aes(fill = rel),
               size = 2, color = 'black',
               alpha = 0.7, shape = 21, stroke = .4) + 
    theme(axis.text.x = element_text(angle = 90, size = 8), 
          legend.position = 'top', panel.background = element_blank(), 
          panel.grid = element_line(colour = 'lightgrey')) +
    scale_fill_viridis_c(option = 'turbo') + 
    labs(x = names(prep.out$control)[2], y = 'BLUP', fill = "Reliability") + 
    guides(fill = guide_colourbar(title.position = 'top', title.hjust = .5, 
                                  barwidth = 9, frame.colour = 'black'))
    
  output$plots$main$DGE.IGE = plott
  
  ### TGV 
  temp = data.frame(gen = main[,1], tgv = main$TGV)
 plott = ggplot(data = temp) + 
    geom_segment(aes(x = stats::reorder(.data$gen, .data$tgv),
                     xend = .data$gen, y = 0, yend = .data$tgv),
                 linewidth = 1) + 
    geom_point(aes(x = stats::reorder(.data$gen, .data$tgv), y = .data$tgv),
               size = 2, color = 'black', alpha = 0.7, shape = 21, stroke = .4, 
               fill = 'darkred') + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, size = 8)) + 
    labs(x = names(prep.out$control)[2], y = 'TGV')
 
 output$plots$main$TGV = plott
  
  if(prep.out$control[,6] > 0) ## Per year ----------------
  {
    output$plots$within = list()
    
    temp = split(within, within[,2])
    
    ## Density
    den.ige.int = do.call(
      ggpubr::ggarrange, 
      c(lapply(temp, function(x){
        aa = ggplot(data = x) + 
          geom_density(aes(x = IGE, after_stat(density), fill = 'Homeostatic'), 
                       alpha = .8)
        d = ggplot2::ggplot_build(aa)$data[[1]]
        aa = aa + geom_area(data = subset(d, x > mean(IGE$IGE) + sd(IGE$IGE)), 
                            aes(x = .data$x, y = .data$y, fill = 'Sensitive'), alpha = 0.8) + 
          geom_area(data = subset(d, x < mean(IGE$IGE) - sd(IGE$IGE)), 
                    aes(x = .data$x, y = .data$y, fill = 'Aggressive'), alpha = .8) +
          geom_density(data = x, aes(x = .data$IGE, after_stat(density)), 
                       color = 'black', linewidth = 1.2, show.legend = F) +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
          labs(x = "IGE", y = 'Density', fill = 'Class', 
               subtitle = paste(names(prep.out$control)[6], unique(x[,2]))) + 
          theme_minimal() +
          theme(legend.position = 'top')
      }), common.legend = TRUE)
    )
    
    output$plots$within$IGE.density = den.ige.int
    
    ### DGE vs IGE
    plott = do.call(
      ggpubr::ggarrange, 
      c(lapply(temp, function(x){
        aa =  ggplot(data = x, aes(x = .data$IGE, y = .data$DGE)) + 
          geom_point(aes(fill = .data$class), colour = 'black', pch = 21, size = 2.2) +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
          labs(fill = 'Class', subtitle = paste(names(prep.out$control)[6],
                                                unique(x[,2]))) +
          theme_minimal() +
          theme(legend.position = 'top') 
      }), common.legend = TRUE)
    )
    
    output$plots$within$DGEvsIGE = plott
    
    ### DGE and IGE
    plott = do.call(
      ggpubr::ggarrange, 
      c(lapply(temp, function(x){
        
        temp2 = data.frame(
          gen = rep(x[, 1], times = 2), 
          blup = c(x$DGE, x$IGE),
          se = c(x$se.DGE, x$se.IGE),
          rel = c(x$rel.DGE, x$rel.IGE),
          comp = rep(c('DGE', 'IGE'), each = prep.out$control[,2])
        )
        
        temp2$gen = factor(temp2$gen, levels = temp2[order(temp2$blup, decreasing = TRUE)[
          order(temp2$blup, decreasing = TRUE) %in% which(temp2$comp == 'DGE')], 1
        ])
        
        aa = ggplot(data = temp2, aes(x = .data$gen, y = .data$blup)) + 
          facet_wrap(.~.data$comp, ncol = 1, scales = 'free_y') + 
          geom_segment(aes(x = .data$gen,
                           xend = .data$gen, y = 0, yend = .data$blup),
                       linewidth = 1) + 
          geom_point(aes(fill = rel),
                     size = 2, color = 'black',
                     alpha = 0.7, shape = 21, stroke = .4) + 
          theme(axis.text.x = element_text(angle = 90, size = 8), 
                legend.position = 'top', panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey')) +
          scale_fill_viridis_c(option = 'turbo') + 
          labs(x = names(prep.out$control)[2], y = 'BLUP', fill = "Reliability") + 
          guides(fill = guide_colourbar(title.position = 'top', title.hjust = .5, 
                                        barwidth = 9, frame.colour = 'black'))
        
      }), common.legend = TRUE)
    )
    
    output$plots$within$DGE.IGE = plott
    
    ### TGV 
    plott = do.call(
      ggpubr::ggarrange, 
      c(lapply(temp, function(x){
        
        aa = data.frame(gen = x[,1], tgv = x$TGV)
        
        bb = ggplot(data = aa) + 
          geom_segment(aes(x = stats::reorder(.data$gen, .data$tgv),
                           xend = .data$gen, y = 0, yend = .data$tgv),
                       linewidth = 1) + 
          geom_point(aes(x = stats::reorder(.data$gen, .data$tgv), y = .data$tgv),
                     size = 2, color = 'black', alpha = 0.7, shape = 21, stroke = .4, 
                     fill = 'darkred') + 
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, size = 8)) + 
          labs(x = names(prep.out$control)[2], y = 'TGV',
               subtitle = paste(names(prep.out$control)[6],
                                unique(x[,2])))
      }))
    )
    
    output$plots$within$TGV = plott

    ## Number of neighbours
    temp = split(prep.out$data, prep.out$data[,names(prep.out$control)[6]])
    
    nneigh = do.call(
      ggpubr::ggarrange,
      c(lapply(temp, function(x){
        Z = x[, 1:prep.out$control[,2]]
        Z = ifelse(crossprod(ifelse(Z == 0, 0, 1)) == 0, 0, 1)
        diag(Z) = 0
        Z = Z[order(rownames(Z)), order(colnames(Z))]
        
        neigh = do.call(rbind, lapply(
          apply(Z, 1, function(y){
            data.frame(
              neigh = names(which(y == 1)),
              class = within[
                which(within[,2] == unique(x[, names(prep.out$control)[6]])),
                'class'
              ][
                match(names(which(y == 1)),
                      within[which(within[,2] == unique(x[, names(prep.out$control)[6]])),
                             c(names(prep.out$control)[2])])
              ]
            )
          }), function(w){
            data.frame(table(w$class))
          }
        ))
        
        neigh$gen = substr(rownames(neigh), start = 1, stop = nchar(unique(rownames(Z))))
        
        ggplot(data = neigh, aes(x = .data$gen, y = .data$Freq, fill = .data$Var1)) +
          geom_bar(stat = "identity", color = 'black') +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, size = 8),
                axis.title = element_text(face = 'bold'),
                legend.position = 'bottom') +
          labs(x = names(prep.out$control)[2], 
               y = "No of different genotypes as neighbours", fill = 'Class', 
               subtitle = paste(names(prep.out$control)[6], 
                                unique(x[,names(prep.out$control)[6]]))) +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6'))
      }), common.legend = TRUE)
    )
    
    output$plots$within$n.neigh = nneigh
    
    ### Grids 
    if(prep.out$control[,7] > 0){
      dat = as.data.frame(model$mf)[,names(prep.out$control)[c(2,4,5,7,6,1)]]
      dat$e = c(model$residuals)
      colnames(dat) = c('gen', 'row', 'col', 'area', 'age', 'y', 'e')
      temp = merge(dat, within[,c(1,2,3,6,9)], by.x = c('gen', 'age'), 
                   by.y = c(names(prep.out$control)[2],
                            names(prep.out$control)[6]))
      
      facet.label.col =  paste(names(prep.out$control)[7],
                           1:as.numeric(prep.out$control[7]))
      names(facet.label.col) = unique(temp$area)
      facet.label.row =  paste(names(prep.out$control)[6],
                               1:as.numeric(prep.out$control[6]))
      names(facet.label.row) = unique(temp$age)
      
      output$plots$within$grid.pheno = ggplot(data = temp, 
                                            aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                fill = .data$y)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), cols = vars(.data$area), 
                   scale = 'free_x', 
                   labeller = labeller(.cols = facet.label.col,
                                       .rows = facet.label.row)) +
        scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'Phenotype') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
      output$plots$within$grid.DGE = ggplot(data = temp, 
                                             aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                 fill = .data$DGE)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), cols = vars(.data$area), 
                   scale = 'free_x', 
                   labeller = labeller(.cols = facet.label.col,
                                       .rows = facet.label.row)) +
        scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'DGE') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
      output$plots$within$grid.IGE = ggplot(data = temp, 
                                            aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                fill = .data$IGE)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), cols = vars(.data$area), 
                   scale = 'free_x', 
                   labeller = labeller(.cols = facet.label.col,
                                       .rows = facet.label.row)) +
        scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'IGE') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
      output$plots$within$grid.IGE.class = ggplot(data = temp, 
                                                  aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                      fill = .data$class)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), cols = vars(.data$area), 
                   scale = 'free_x', 
                   labeller = labeller(.cols = facet.label.col,
                                       .rows = facet.label.row)) +
        scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'IGE: Classes') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'bottom', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
    }else{
      dat = as.data.frame(model$mf)[,names(prep.out$control)[c(2,4,5,6,1)]]
      dat$e = c(model$residuals)
      colnames(dat) = c('gen', 'row', 'col', 'age', 'y', 'e')
      temp = merge(dat, within[,c(1,2,3,6,9)], by.x = c('gen', 'age'), 
                   by.y = c(names(prep.out$control)[2],
                            names(prep.out$control)[6]))
      
      facet.label.row =  paste(names(prep.out$control)[6],
                               1:as.numeric(prep.out$control[6]))
      names(facet.label.row) = unique(temp$age)
      
      output$plots$within$grid.pheno = ggplot(data = temp, 
                                              aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                  fill = .data$y)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), scale = 'free_x', 
                   labeller = labeller(.rows = facet.label.row)) +
        scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'Phenotype') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
      output$plots$within$grid.DGE = ggplot(data = temp, 
                                            aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                fill = .data$DGE)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), scale = 'free_x', 
                   labeller = labeller(.rows = facet.label.row)) +
        scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'DGE') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
      output$plots$within$grid.IGE = ggplot(data = temp, 
                                            aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                fill = .data$IGE)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), scale = 'free_x', 
                   labeller = labeller(.rows = facet.label.row)) +
        scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'IGE') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
      
      output$plots$within$grid.IGE.class = ggplot(data = temp, 
                                                  aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                                                      fill = .data$class)) + 
        geom_tile(color = 'black') + 
        facet_grid(rows = vars(.data$age), scale = 'free_x', 
                   labeller = labeller(.rows = facet.label.row)) +
        scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
        labs(x = "Row", y = 'Column', fill = '',
             title = 'IGE: Classes') + 
        theme(plot.title = element_text(hjust = .5), legend.position = 'bottom', 
              panel.background = element_blank(), 
              panel.grid = element_line(colour = 'lightgrey'))
    }
    
  }else{
    ## Number of neighbours
    Z = prep.out$Z
    Z = ifelse(crossprod(ifelse(Z == 0, 0, 1)) == 0, 0, 1)
    diag(Z) = 0
    Z = Z[order(rownames(Z)), order(colnames(Z))]
    
    neigh = do.call(rbind, lapply(
      apply(Z, 1, function(y){
        data.frame(
          neigh = names(which(y == 1)),
          class = main[match(names(which(y == 1)), main[,1]), 8]
        )
      }), function(w){
        data.frame(table(w$class))
      }
    ))
    
    neigh$gen = substr(rownames(neigh), start = 1, stop = nchar(unique(rownames(Z))))
    
    nneigh = ggplot(data = neigh, aes(x = .data$gen, y = .data$Freq, 
                                      fill = .data$Var1)) +
      geom_bar(stat = "identity", color = 'black') +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, size = 8),
            axis.title = element_text(face = 'bold'),
            legend.position = 'bottom') +
      labs(x = names(prep.out$control)[2], 
           y = "Number of neighbours", fill = 'Class') +
      scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6'));nneigh
    
    output$plots$main$n.neigh = nneigh
    
    
    ### Grids 
    if(prep.out$control[,7] > 0){
      
      dat = as.data.frame(model$mf)[,names(prep.out$control)[c(2,4,5,7,1)]]
      dat$e = c(model$residuals)
      colnames(dat) = c('gen', 'row', 'col', 'area', 'y', 'e')
      temp = merge(dat, main[,c(1,2,5,8)], by.x = 'gen', by.y = names(prep.out$control)[2])
      
      facet.label =  paste(names(prep.out$control)[7],
                           1:as.numeric(prep.out$control[7]))
      names(facet.label) = unique(temp$area)
      
      plott = ggpubr::ggarrange(
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$y)) + 
          geom_tile(color = 'black') + 
          facet_wrap(.~.data$area, scale = 'free_x', 
                     labeller = labeller(.cols = facet.label)) +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = '',
               title = 'Phenotype') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey')),
        # ggplot(data = temp, 
        #        aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
        #            fill = .data$e)) + 
        #   geom_tile(color = 'black') + 
        #   facet_wrap(.~.data$area, scale = 'free_x', 
        #              labeller = labeller(.cols = facet.label)) +
        #   scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        #   labs(x = "Row", y = 'Column', fill = '',
        #        title = 'Residue') + 
        #   theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
        #         panel.background = element_blank(), 
        #         panel.grid = element_line(colour = 'lightgrey')), 
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$DGE)) + 
          geom_tile(color = 'black') + 
          facet_wrap(.~.data$area, scale = 'free_x', 
                     labeller = labeller(.cols = facet.label)) +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = '',
               title = 'DGE') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey')),
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$class)) + 
          geom_tile(color = 'black') + 
          facet_wrap(.~.data$area, scale = 'free_x', 
                     labeller = labeller(.cols = facet.label)) +
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
          labs(x = "Row", y = 'Column', fill = '',
               title = 'IGE: Classes') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'top', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey')),
        ggplot(data = temp, 
               aes(x = as.numeric(.data$row), y = as.numeric(.data$col), 
                   fill = .data$IGE)) + 
          geom_tile(color = 'black') + 
          facet_wrap(.~.data$area, scale = 'free_x', 
                     labeller = labeller(.cols = facet.label)) +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = '',
               title = 'IGE') + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right', 
                panel.background = element_blank(), 
                panel.grid = element_line(colour = 'lightgrey'))
      )
      
      output$plots$main$grid = plott
      
    }else{
      temp = merge(
        data.frame(prep.out$neigh_check[,1:4], e = model$residuals),
        main[,c(1,2,5,8)], by.x = 'trat', by.y = names(prep.out$control)[2]
      )
      
      plott = ggpubr::ggarrange(
        ggplot(data = temp, 
               aes(x = .data$row, y = .data$col, fill = .data$y_focal)) + 
          geom_tile(color = 'black') +
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = '',
               title = 'Phenotype') + 
          theme_minimal() + 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right'),
        # ggplot(data = temp, 
        #        aes(x = .data$row, y = .data$col, fill = .data$e)) + 
        #   geom_tile(color = 'black') + 
        #   scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
        #   labs(x = "Row", y = 'Column', fill = '', title = 'Residual') + 
        #   theme_minimal()+ 
        #   theme(plot.title = element_text(hjust = .5)), 
        ggplot(data = temp, 
               aes(x = .data$row, y = .data$col, fill = .data$DGE)) + 
          geom_tile(color = 'black') + 
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = '', title = 'DGE') + 
          theme_minimal()+ 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right'),
        ggplot(data = temp, 
               aes(x = .data$row, y = .data$col, fill = .data$class)) + 
          geom_tile(color = 'black')+
          scale_fill_manual(values = c('#fd8d3c', '#edf8e9', '#6baed6')) + 
          labs(x = "Row", y = 'Column', fill = '', title = 'IGE: Classes') + 
          theme_minimal()+ 
          theme(plot.title = element_text(hjust = .5), legend.position = 'top'),
        ggplot(data = temp, 
               aes(x = .data$row, y = .data$col, fill = .data$IGE)) + 
          geom_tile(color = 'black') + 
          scale_fill_viridis_c(option = 'turbo', na.value = 'white') + 
          labs(x = "Row", y = 'Column', fill = '', title = 'IGE') + 
          theme_minimal()+ 
          theme(plot.title = element_text(hjust = .5), legend.position = 'right')
      )
      
      output$plots$main$grid = plott
    }
  }
  remove(prep.out, envir = .GlobalEnv)
  
  return(output)
}









