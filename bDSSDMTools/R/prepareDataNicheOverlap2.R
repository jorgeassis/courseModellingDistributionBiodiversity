#' Prepare ecospat Niche Grids for Two Species
#'
#' Samples presence cells from two binary rasters, extracts the environment at
#' those cells, and summarises the two niches on the first two axes of a PCA.
#' The PCA is fitted to the pooled sample of both species. The same scores are
#' then used as the species occurrences and as the background passed to
#' `ecospat::ecospat.grid.clim.dyn`.
#'
#' Up to `maxPoints` cells are drawn at random from each species, without
#' replacement. Set the random seed before calling this function if the sample
#' must be repeated. A contribution plot of the PCA axes is drawn.
#'
#' @param environment A terra `SpatRaster` with at least two environmental layers.
#' @param rasterSPA A single-layer terra `SpatRaster`. Cells with value `1` are
#'   presences of the first species.
#' @param rasterSPB A single-layer terra `SpatRaster` for the second species.
#' @param maxPoints Maximum number of presence cells kept for each species.
#'   Defaults to `1000`.
#'
#' @return A named list with the dynamic occurrence grids of species A
#'   (`grid.clim.sp.1`, from `rasterSPA`) and species B (`grid.clim.sp.2`, from
#'   `rasterSPB`). Each element is the object returned by
#'   `ecospat::ecospat.grid.clim.dyn`.
#' @export
#' @importFrom terra as.data.frame extract nlyr
#' @importFrom stats complete.cases var
#' @importFrom ade4 dudi.pca suprow
#' @importFrom ecospat ecospat.plot.contrib ecospat.grid.clim.dyn
#' @importFrom graphics par
#' @seealso [prepareDataNicheOverlap()]
#'
#' @examples
#' \dontrun{
#' grids <- prepareDataNicheOverlap2(environment, rasterSPA, rasterSPB)
#' names(grids)
#' }
prepareDataNicheOverlap2 <- function(environment, rasterSPA, rasterSPB, maxPoints = 1000) {

  if (!is.numeric(maxPoints) || length(maxPoints) != 1 || maxPoints < 2) {
    stop("Error :: maxPoints must be a number greater than 1")
  }
  if (terra::nlyr(environment) < 2) {
    stop("Error :: environment must have at least two layers")
  }

  environment.used.by.species.1 <- .sampleRows(
    .environmentAtPresences(environment, rasterSPA),
    maxPoints
  )
  environment.used.by.species.2 <- .sampleRows(
    .environmentAtPresences(environment, rasterSPB),
    maxPoints
  )
  if (nrow(environment.used.by.species.1) <= terra::nlyr(environment) ||
      nrow(environment.used.by.species.2) <= terra::nlyr(environment)) {
    stop("Error :: after sampling, each species still needs more presence cells than environmental layers")
  }

  env.used <- rbind(environment.used.by.species.1, environment.used.by.species.2)
  constant <- vapply(env.used, function(z) {
    v <- stats::var(z)
    !is.finite(v) || v == 0
  }, logical(1))
  if (any(constant)) {
    stop(
      "Error :: environmental layers with no variance at the sampled presences: ",
      paste(names(env.used)[constant], collapse = ", ")
    )
  }

  pca.env.used <- ade4::dudi.pca(env.used, scannf = FALSE, nf = 2)
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  ecospat::ecospat.plot.contrib(contrib = pca.env.used$co, eigen = pca.env.used$eig)

  scores.globclim <- pca.env.used$li
  scores.sp.1 <- ade4::suprow(pca.env.used, Xsup = environment.used.by.species.1)$li
  scores.sp.2 <- ade4::suprow(pca.env.used, Xsup = environment.used.by.species.2)$li

  # The session workflow uses the presence environments both as the species
  # records and as the background climate of each species.
  grid.clim.sp.1 <- ecospat::ecospat.grid.clim.dyn(
    glob = scores.globclim,
    glob1 = scores.sp.1,
    sp = scores.sp.1,
    R = 100,
    th.sp = 0
  )
  grid.clim.sp.2 <- ecospat::ecospat.grid.clim.dyn(
    glob = scores.globclim,
    glob1 = scores.sp.2,
    sp = scores.sp.2,
    R = 100,
    th.sp = 0
  )

  list(grid.clim.sp.1 = grid.clim.sp.1, grid.clim.sp.2 = grid.clim.sp.2)

}
