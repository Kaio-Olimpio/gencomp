% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/composite.R
\name{composite}
\alias{composite}
\title{Simulates a clonal composite}
\usage{
composite(
  prep.out,
  model,
  resp.out,
  d.row.col,
  d.weight = TRUE,
  selected,
  nsim = 10,
  verbose = TRUE
)
}
\arguments{
\item{prep.out}{A \code{comprepfor} object.}

\item{model}{A \code{compmod} object obtained using \code{\link[=asr]{asr()}} or \code{\link[=asr_ma]{asr_ma()}}.}

\item{resp.out}{A \code{comresp} object.}

\item{d.row.col}{A vector of size two. The first element contain the distance between
rows, and second the distance between columns of the simulated grid.}

\item{d.weight}{A logical value. If \code{TRUE} (default) the predicted mean
of each plant in the grid will be weighted by the inverse of the distance between rows,
columns and diagonals.}

\item{selected}{A vector with the names of the clones selected to compose the clonal composite.
The names must be present in the \code{comresp} object.}

\item{nsim}{An integer defining the number grid simulations. If \code{nsim > 1},
the function will estimate the 95\% confidence interval of the predicted means using
a bootstrap process. Defaults to 10.}

\item{verbose}{A logical value. If \code{TRUE}, shows the function progress. Defaults to \code{FALSE}.}
}
\value{
The function returns a data.frame with the predicted mean of each clone
and their respective 95\% confidence interval. The composite performance is obtained by
averaging the predicted mean of all clones.
}
\description{
This function simulates clonal composites using outputs of a genetic-spatial
competition model fitted using \code{\link[=asr]{asr()}} or \code{\link[=asr_ma]{asr_ma()}}. Currently,
only forestry data is accepted.
}
\details{
Considering the direct (DGE) and indirect genetic effects (IGE) of the selected clones,
the function simulates grids. Clones are positioned differently in each simulation,
which enables the modification of focal tree-neighbour dynamics. In each simulation,
the expected mean of each clone is predicted using the following equation (Ferreira et al., 2023):
\deqn{\hat{y}_{ij} = \hat{\mu} + \hat{d}_i + \sum^n_{i \neq j}{\hat{c}_j}}
where \eqn{d_i} is the DGE of the i\if{html}{\out{<sup>}}th\if{html}{\out{</sup>}} focal tree, and
\eqn{c_j} is the IGE of the j\if{html}{\out{<sup>}}th\if{html}{\out{</sup>}} neighbour. If \code{d.weight = TRUE}, the
IGE is divided by the distance between the focal tree and its neighbours:
\deqn{\hat{y}_{ij} = \hat{\mu} + \hat{d}_i + \sum^n_{i \neq j}{\frac{1}{dist_{ij}} \times \hat{c}_j}}
}
\examples{
\donttest{
library(gencomp)
 comp_mat = prepfor(data = euca, gen = 'clone', area = 'area',
                   plt = 'tree', age = 'age', row = 'row', col = 'col',
                   dist.col = 3, dist.row = 2, trait = 'MAI', method = 'SK',
                   n.dec = 3, verbose = FALSE, effs = c("block"))
 model = asr_ma(prep.out = comp_mat,
                fixed = MAI~ age, 
                random = ~ block:age, 
                lrtest = TRUE, 
                spatial = TRUE, 
                cor = TRUE, 
                maxit = 20)
             
 results = resp(prep.out = comp_mat, model = model, weight.tgv = FALSE, sd.class = 1)
 
 cc = composite(prep.out = comp_mat, model = model, resp.out = results,
                d.row.col = c(3,3), d.weight = TRUE, nsim = 10, verbose = TRUE, 
                selected = results$blups$main[order(results$blups$main$TGV, 
                                                    decreasing = TRUE),1][1:10])
 }
}
\references{
Ferreira, F.M., Chaves, S.F.S., Bhering, L.L., Alves, R.S., Takahashi, E.K.,
Sousa, J.E., Resende, M.D.V., Leite, F.P., Gezan, S.A., Viana, J.M.,
Fernandes, S.B., Dias, K.O.G., 2023. A novel strategy to predict clonal composites
by jointly modeling spatial variation and genetic competition. Forest Ecology
and Management 548, 121393. \doi{https://doi.org/10.1016/j.foreco.2023.121393}
}
\seealso{
\link{prepfor}, \link{asr}, \link{asr_ma}, \link{resp}
}
