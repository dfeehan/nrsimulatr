###############################################################################
#' create a social network from an existing igraph object
#' 
#' Given an existing igraph object, add attributes for the frame
#' population and for the hidden population.
#'
#'
#' @param params an optional list containing hidden.pop and frame.pop
#' @return the \code{igraph} object with each vertex's hidden and frame pop membership 
#' indicated by in.F and in.H
#' @export
generate_graph.existing_igraph <- function(params, g) {
  
  ## these are the names of vertex attributes (assumed already to exist)
  ## that will be reported
  groups <- params$groups
  
  ## the names of the groups that are in the frame population
  gps.in.F <- params$gps.in.F
  
  ## the names of the groups that are in the hidden population
  gps.in.H <- params$gps.in.H
  
  # label vertices w/ dummy variables for membership in F
  V(g)$in.F <- 0
  
  for (gp in gps.in.F) {
     g <- set_vertex_attr(g, 
                          'in.F', 
                          ## NB: we're assuming that the attributes are TRUE/FALSE or 1/0
                          index=as.logical(vertex_attr(g, gp)), 
                          value=1) 
  }
  
  # ... and for membership in H
  V(g)$in.H <- 0
  
  for (gp in gps.in.H) {
    g <- set_vertex_attr(g, 
                         'in.H', 
                         ## NB: we're assuming that the attributes are TRUE/FALSE or 1/0
                         index=as.logical(vertex_attr(g, gp)), 
                         value=1) 
  }
  
  ## calculate the degree and add it as an attribute
  ## note that this is subtly different from edge-reports, 
  ## where y.degree and v.degree can be affected
  ## by reporting errors
  g <- set.vertex.attribute(g, 
                            'degree', 
                            value=degree(g))
  
  g <- set.graph.attribute(g, "groups", groups)
  g <- set.graph.attribute(g, "gps.in.F", gps.in.F)
  g <- set.graph.attribute(g, "gps.in.H", gps.in.H)
  
  class(g) <- append(class(g), "existing_igraph")
  
  return(g)
  
}

###############################################################################
#' create a parameter object for an existing igraph
#' 
#' Given an existing igraph object, add attributes for the frame
#' population and for the hidden population.
#' 
#' The \code{params} vector should have at least three entries:
#' \itemize{
#' \item \code{groups} is a vector that has a list of attributes that
#' each node in the \code{igraph} object \code{g} has. These attributes
#' should all have values of 0/1 or TRUE/FALSE.
#' \item \code{gps.in.F} is a vector that has a list of group names; nodes in
#' those groups are taken to be in the frame population
#' \item \code{gps.in.H} is a vector that has a list of group names; nodes in
#' those groups are taken to be in the hidden population
#' }'
#' 
#' @param params the parameters
#' @return an object with the parameters
#' @examples 
#' my_params <- existing_params(list(group.names=c('a', 'b', 'c', 'd'),
#'                                   gps.in.F=c('b', 'c'),
#'                                   gps.in.H=c('a', 'b')))
#' @export
existing_params <- function(params) {
  
  # TODO - eventually, check that the params are all there?
  
  class(params) <- "existing_igraph"
  
  return(params)
}


###############################################################################
#' for each vertex, compute the number of  connections to each block
#'
#' If this function is called on an undirected social network graph (with mode='all'),
#' then it counts the true number of social connections between each node and
#' the various blocks.
#'
#' On the other hand, if this function is called on a directed reporting graph,
#' then it computes the reported number of social connections (with mode='in')
#' or visibility (with mode='out').
#'
#' @details
#' 
#' This function is useful for computing edge counts (i.e. )
#' between
#' individual vertices and the blocks. These edge counts are the building
#' blocks of many network reporting estimators.
#' 
#' Since this function is based on the stochastic block model, it assumes
#' that the groups are mutually exclusive. This could be modified in the
#' future.
#'
#' Note that, counter-intuitively, mode="in" will compute out-reports and
#' mode="out" will compute in-reports. This is the \code{igraph} convention.
#' For undirected graphs (for example, when computing degrees in the social
#' network), use "all".
#' 
#' For example, if the prefix is "y." and the mode is "in", we are asking
#' \code{report.sbm.edges} to count each vertex's number of out-reports
#' about each block. If we have three blocks, A, B, and C, then in the
#' graph that \code{report.sbm.edges} returns, each vertex will have
#' four new attributes: \code{y.degree}, \code{y.A}, \code{y.B}, and \code{y.C}.
#'
#' @param g the \code{igraph} object
#' @param prefix the prefix to use in the variable names that are attached
#'        (useful if this function will be used to compute reports more than once)
#' @param mode one of "all", "in", or "out". see Details
#' @return the \code{igraph} object with new attributes affixed to each vertex
#'         (see Details)
#' @export
report_edges.existing_igraph <- function(g, prefix='d.', mode="all") {
  
  # compute the degree of each vertex
  # note that for undirected graphs, mode='out' has no effect; this just
  # returns the degree
  g <- set.vertex.attribute(g, paste0(prefix, 'degree'), value=degree(g, mode=mode))
  
  adj.matrix <- get.adjacency(g)
  
  groups <- get.graph.attribute(g, 'groups')
  
  
  ## go through and get reported connections to each block
  for (gp in groups) {
    
    in.gp <- as.numeric(get.vertex.attribute(g, gp)) 
    
    res.name <- paste0(prefix, gp)
    
    if (mode=='out') {
      # in this case, we want visibilities
      # which are A^T x; we'll compute them via
      # we'll compute x^T A instead
      num.rep.conns <- as.numeric(t(in.gp) %*% adj.matrix)
    } else {
      num.rep.conns <- as.numeric(adj.matrix %*% in.gp)
    }
    
    g <- set.vertex.attribute(g, 
                              res.name,
                              value=num.rep.conns)
  }
  
  return(g)
  
}