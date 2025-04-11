#' List Multiple Bio-Oracle Layers
#'
#' List specified Bio-Oracle layers
#'
list_multiple_layers <- function() {

  # list available layers in Bio-ORACLE
  layersList <- biooracler::list_layers()
  return(layersList)

}
