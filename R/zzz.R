.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "MosaiClusteR ", utils::packageVersion("MosaiClusteR"),
    " - an umbrella framework for multi-source clustering. ",
    "See vignette('MosaiClusteR') and ?M_ABCpp.")
}
