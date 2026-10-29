# Names used inside ggplot2 non-standard evaluation. Declaring them here keeps
# R CMD check from treating those plot columns as missing objects.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Lat", "Lon", "Predictor", "Set", "Variable",
    "f_calibrated_score", "folds", "relImportance", "relImportance.sd",
    "value", "variable", "x", "xa", "y", "ya"
  ))
}
