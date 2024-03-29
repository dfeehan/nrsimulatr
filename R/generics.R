###############################################################################
#' generate a social network
#' 
#' this is a generic function to generate a social network
#'
#' @param params object with the parameters governing how the network
#'        should be generated
#' @return the \code{igraph} object with the generated social network
#' @export
generate_graph <- function(params, ...) {
    UseMethod("generate_graph")
}

###############################################################################
#' make a reporting graph from a social network (ARD)
#' 
#' this is a generic function to generate a reporting graph with
#' aggregated relational data from
#' reporting parameters and a social network
#'
#' @param params object with the parameters governing reporting (i.e., how the 
#'        reporting network should be generated
#' @return the \code{igraph} object with the generated social network
#' @export
reporting_graph <- function(params, g, ...) {
    UseMethod("reporting_graph", params)
}

###############################################################################
#' report edges (ARD)
#' 
#' this is a generic function to report edges as ARD from a social network or
#' a reporting graph
#'
#' @param g the social network or reporting graph
#' @return the \code{igraph} object with the generated social network
#' @export
report_edges <- function(g, ...) {
    UseMethod("report_edges")
}


###############################################################################
#' make a detailed reporting graph from a social network
#' 
#' this is a generic function to generate a reporting graph with
#' detailed reports from
#' reporting parameters and a social network
#'
#' @param params object with the parameters governing reporting (i.e., how the 
#'        reporting network should be generated
#' @return the \code{igraph} object with the generated social network
#' @export
reporting_graph_detailed <- function(params, g, ...) {
  UseMethod("reporting_graph_detailed", params)
}

###############################################################################
#' report detailed edges
#' 
#' this is a generic function to report edges one by one from a social network 
#' or a reporting graph
#'
#' @param g the social network or reporting graph
#' @return the \code{igraph} object with the generated social network
#' @export
report_detailed_edges <- function(g, ...) {
  UseMethod("report_detailed_edges")
}




###############################################################################
#' produce dataset from a sample of a reporting graph or social network
#' 
#' this is a generic function to generate a dataset with information
#' from a sample of nodes in a reporting graph or social network.
#' the parameters describe how the sample is obtained.
#'
#' @param params object with the parameters governing how the sample should
#'        be obtained
#' @param g the social network or reporting graph to sample' 
#' @return dataframe with the attributes of the sampled vertices
#' @export
sample_graph <- function(params, g, ...) {
    UseMethod("sample_graph", params)
}


