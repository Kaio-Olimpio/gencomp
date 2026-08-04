rm(list=ls())

# Library -----------------------------------------------------------------
library(dplyr)
library(devtools)
library(MASS)
library(FieldSimR)
library(asreml)

# Simulation parameters ---------------------------------------------------
s2d = 25
s2c = 5
rdc = -.8
sdc = rdc*sqrt(s2d)*sqrt(s2c)
h2 = .4
s2e = s2d/h2 - s2d
drow = 2
dcol = 3
mu = 25
seed = 1997

# Founders ----------------------------------------------------------------
set.seed(seed * 2)
nfound = 20
C = matrix(c(s2d, rdc*sqrt(s2d)*sqrt(s2c), rdc*sqrt(s2d)*sqrt(s2c), s2c), 
           nrow= 2, byrow = TRUE)
founders = data.frame(
  founder = paste0("F", sprintf(paste0('%0', nchar(nfound),'d'), seq(1:nfound))),
  mvrnorm(n = nfound, mu = c(0,0), Sigma = C)
) |> rename(direct = X1, competition = X2)

diallel = matrix(NA, nrow = nfound, ncol = nfound, 
                    dimnames = list(founders$founder, founders$founder))
for (i in 1:nrow(diallel)) {
  aux = seq(i,i+3)[-1]
  if(any(aux > nfound)) aux[which(aux > nfound)] = aux[which(aux > nfound)] - nfound
  diallel[i, aux] = 1
}
diallel = cbind(rownames(diallel), diallel)
nfam = sum(diallel == "1", na.rm = TRUE)

# Clonal progeny trials ======

## Simulation of direct and indirect effects ----
nprog = 100  ## Simulating 100 progenies per family
simu_gen = do.call(rbind, apply(diallel, 1, function(x){
  d_Gfp = founders[which(founders$founder == x[1]),"direct"]
  c_Gfp = founders[which(founders$founder == x[1]),"competition"]
  d_Gmp = founders[which(founders$founder %in% names(na.exclude(x))[-1]),"direct"]
  c_Gmp = founders[which(founders$founder %in% names(na.exclude(x))[-1]),"competition"]
  
  aux = rbind(
    .5*d_Gfp + .5*d_Gmp,
    .5*c_Gfp + .5*c_Gmp
  )
  colnames(aux) = paste(x[1], names(na.exclude(x))[-1], sep="x")
  rownames(aux) = c("d", "c")
  set.seed(as.numeric(gsub("F", "", x[1])) * seed)
  aux = apply(aux, 2, function(w){
    mvrnorm(nprog, mu = w, Sigma = .5 * C)
  })
  aux = as.data.frame(aux)
  aux$effect = rep(c("d", "c"), each = nprog)
  dd = aux[which(aux$effect == "d"),]
  cc = aux[which(aux$effect == "c"),]
  aux = merge(
    reshape(dd, direction = "long", varying = list(1:3), times = colnames(dd)[-4],
            v.names = "direct"),
    reshape(cc, direction = "long", varying = list(1:3), times = colnames(cc)[-4],
            v.names = "competition"),
    by = c("time", "id")
  )
  aux$id = as.numeric(aux$id)
  aux = aux[,c(1,2,4,6)]
  colnames(aux)[1] = "prog"
  aux = aux[order(aux$prog, aux$id),]
  return(aux)
}))

## Simulation of residuals ----
#' I will sample 10 trees per family, then repeat each of them five times
#' 60 families x 10 trees x 5 reps = 3000 plots
nplots = 3000
# do.call(rbind, findpairs(nplots))

nrow = 40
ncol = 75
numrept = 5
prop_spatial = 0.3 # proportion of spatial error variance
prop_ext = 0.01 # proportion of extraneous error variance
1 - (prop_spatial + prop_ext) # proportion of random error variance

set.seed(seed * 6)
error_df1 = field_trial_error(
  ntraits = 1,
  nenvs = 1,
  nblocks = numrept,
  block.dir = "row",
  ncols = ncol,
  nrows = nrow,
  varR = s2e,
  spatial.model = "AR1",
  col.cor = -.25,
  row.cor = -.5,
  prop.spatial = prop_spatial,
  ext.ord = "zig-zag",
  ext.dir = "col",
  prop.ext = prop_ext,
  return.effects = TRUE
)

set.seed(seed * 7)
temp = do.call(rbind, lapply(split(simu_gen, simu_gen$prog), function(x) {
  x[sample(1:nrow(x), 10), ]
}))
simu_df = do.call(rbind, 
        lapply(split(error_df1$error.df, error_df1$error.df$block), function(x){
          cbind(x, temp[sample(1:nrow(temp)),])
        }))
simu_df = simu_df[,c("prog", "id", "col", "row", "block", "direct", "competition", "e.Trait1")]
simu_df$row = as.numeric(simu_df$row)
simu_df$col = as.numeric(simu_df$col)
simu_df = simu_df[order(simu_df$row, simu_df$col), ]
simu_df$id = paste(simu_df$prog, simu_df$id, sep = "_")
rownames(simu_df) = NULL

## Simulation of phenotypes (using the competition matrix) -----
N = nrow(simu_df)
id_levels = unique(simu_df$id)
n_levels = length(id_levels)

Zs = matrix(0, nrow = N, ncol = n_levels, dimnames = list(1:N, id_levels))

p = dcol / drow
pheno = numeric(N)

for (i in 1:N) {
  r_i = simu_df$row[i]
  c_i = simu_df$col[i]

  idx_diag = which(simu_df$row %in% (r_i + c(-1, 1)) & simu_df$col %in% (c_i + c(-1, 1)))
  idx_col  = which(simu_df$col == c_i & simu_df$row %in% (r_i + c(-1, 1)))
  idx_row  = which(simu_df$row == r_i & simu_df$col %in% (c_i + c(-1, 1)))
  
  n_d = length(idx_diag)
  n_c = length(idx_col)
  n_r = length(idx_row)
  
  if (n_d + n_c + n_r > 0) {
    f_D = round(p / sqrt((n_r * p^4) + (n_r * p^2) + (n_c * p^2) + (n_d * p^2) + n_c), 4)
    f_C = round((f_D * sqrt(1 + p^2)) / p, 4)
    f_R = round(f_D * sqrt(1 + p^2), 4)
    
    if (n_d > 0) {
      for (k in idx_diag) {
        g = as.character(simu_df$id[k])
        Zs[i, g] = Zs[i, g] + f_D
      }
    }
    
    if (n_c > 0) {
      for (k in idx_col) {
        g = as.character(simu_df$id[k])
        Zs[i, g] = Zs[i, g] + f_C
      }
    }
    
    if (n_r > 0) {
      for (k in idx_row) {
        g = as.character(simu_df$id[k])
        Zs[i, g] = Zs[i, g] + f_R
      }
    }
  }
    comp_eff = sum(Zs[i, ] * simu_df$competition[match(id_levels, simu_df$id)])
    pheno[i] = mu + simu_df$direct[i] + comp_eff + simu_df$e.Trait1[i]
}
simu_df$pheno = pheno

ped = unique(as.data.frame(cbind(
  do.call(rbind, strsplit(as.character(simu_df[, 1]), split = "x")), as.character(simu_df[, 2])
)))
ped = ped[,c(3,1,2)]
ainv = ainverse(ped)

simu_df = transform(
  simu_df,
  id = as.factor(id),
  prog = as.factor(prog),
  col = as.factor(col),
  row = as.factor(row),
  block = as.factor(block)
)

simu_df = cbind(matrix(
  0,
  ncol = 20,
  nrow = nrow(simu_df),
  dimnames = list(NULL, paste0("F", sprintf(
    paste0('%0', nchar(nfound), 'd'), seq(1:nfound)
  )))
), Zs, simu_df)
trat = attr(ainv, 'rowNames')

colnames(simu_df)[1:length(trat)]
simu_df[,1:length(trat)] = simu_df[,1:length(trat)][,match(trat, colnames(simu_df)[1:length(trat)])]

## Models -----
mod1.1 = asreml(fixed = pheno ~ block, 
              random = ~ id, 
              data = simu_df)
summary(mod1.1)$varcomp

mod1.2 = asreml(fixed = pheno ~ block, 
                random = ~ vm(id, ainv), 
                data = simu_df)
summary(mod1.2)$varcomp



mod2.1 = asreml(fixed = pheno ~ block, 
              random = ~ id,
              residual = ~ar1v(row):ar1(col), 
              data = simu_df)
summary(mod2.1)$varcomp

mod2.2 = asreml(fixed = pheno ~ block, 
              random = ~  vm(id, ainv),
              residual = ~ar1v(row):ar1(col), 
              data = simu_df)
summary(mod2.2)$varcomp


mod3.1 = asreml(fixed = pheno ~ block, 
                random = ~ str(~id + grp(g1), ~corh(2):id(id)),
                residual = ~ar1v(row):ar1(col), 
                data = simu_df,
                group = list(g1 = 21:620),
                maxit = 40)
summary(mod3.1)$varcomp

mod3.2 = asreml(fixed = pheno ~ block, 
                random = ~ str(~vm(id, ainv) + grp(g1), ~corh(2):vm(id, ainv)),
                residual = ~ar1v(row):ar1(col), 
                data = simu_df,
                group = list(g1 = 1:620),
                maxit = 40)
summary(mod3.2)$varcomp


lrt(
  asreml(fixed = pheno ~ block, 
         random = ~ id + grp(g1),
         residual = ~ar1v(col):ar1(row), 
         data = simu_df,
         group = list(g1 = 21:620),
         maxit = 40),
  asreml(fixed = pheno ~ block, 
         random = ~ id,
         residual = ~ar1v(col):ar1(row), 
         data = simu_df,
         group = list(g1 = 21:620),
         maxit = 40)
)

lrt(
  asreml(fixed = pheno ~ block, 
         random = ~ vm(id, ainv) + grp(g1),
         residual = ~ar1v(col):ar1(row), 
         data = simu_df,
         group = list(g1 = 1:620),
         maxit = 40),
  asreml(fixed = pheno ~ block, 
         random = ~ vm(id, ainv),
         residual = ~ar1v(col):ar1(row), 
         data = simu_df,
         group = list(g1 = 1:620),
         maxit = 40)
)

# BLUPs -------------------------------------------------------------------
blu = summary(mod3.2, coef = TRUE)$coef.random
blup = data.frame(d_blup = blu[grep("vm", rownames(blu)), 1]) |> rownames_to_column("id") |>
  mutate(id = gsub("vm\\(id, ainv\\)_", "", id)) |> 
  left_join(data.frame(c_blup = blu[grep("grp", rownames(blu)), 1]) |> rownames_to_column("id") |>
              mutate(id = gsub("grp\\(g1\\)_", "", id))) |> 
  left_join(unique(simu_df[,c("id", "direct", "competition")])) |> 
  left_join(founders, by = c("id" = "founder")) |> 
  mutate(direct = ifelse(is.na(direct.x), direct.y, direct.x),
         competition = ifelse(is.na(competition.x), competition.y, competition.x)) |> 
  dplyr::select(-direct.x,-direct.y,-competition.x,-competition.y)


ggplot(data = blup, aes(x = d_blup, y = direct)) + 
  geom_point(alpha = .5, size = 2) + 
  theme_bw()

ggplot(data = blup, aes(x = c_blup, y = competition)) + 
  geom_point(alpha = .5, size = 2) + 
  theme_bw()

ggplot(data = blup, aes(x = d_blup, y = c_blup)) + 
  geom_point(alpha = .5, size = 2) + 
  theme_bw()

# Testing gencomp ---------------------------------------------------------
rm(list=ls())
devtools::load_all()
devtools::document()

comp = prepfor(
  data = cpt,
  gen = "id",
  row = "row",
  col = "col",
  trait = "pheno",
  effs = "block",
  dist.row = 2,
  dist.col = 3,
  verbose = TRUE,
  n.dec = 4
)

ped = unique(as.data.frame(cbind(
  do.call(rbind, strsplit(as.character(cpt[, 1]), split = "x")), as.character(cpt[, 2])
)))
ped = ped[,c(3,1,2)]
ainv = ainverse(ped)

mod = asr(
  prep.out = comp,
  fixed = pheno ~ block,
  random = ~ 1,
  spatial = TRUE,
  cor = TRUE,
  lrtest = TRUE,
  K = ainv
)
summary(mod)$varcomp
mod$lrt

result = resp(
  prep.out = comp,
  model = mod,
  weight.tgv = TRUE,
  sd.class = 1
)

devtools::check()


