###############################################################################
#' generate a social network
#' 
#' this is a generic function to generate a social network
#'
#' @param params object with the parameters governing how the network
#'        should be generated
#' @return the \code{igraph} object with the generated social network
#' @export
generate_graph <- function(params) {
    UseMethod("generate_graph")
}

reporting_graph <- function(params, g, ...) {
    UseMethod("reporting_graph", params)
}

report_edges <- function(g, ...) {
    UseMethod("report_edges")
}

sample_graph <- function(params, g, ...) {
    UseMethod("sample_graph", params)
}


