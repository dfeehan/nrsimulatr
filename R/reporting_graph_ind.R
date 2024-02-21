###############################################################################
# functions related to individual reporting
###############################################################################


###############################################################################
#' Create a set of parameters describing reporting behavior that is not perfectly accurate
#' and that varies by individual
#'
#'
#' @details
#' For now, the only simulation parameter that makes a difference is
#' \eqn{\tau}.
#'
#' @param params a list of params, which must contain tau
#' @return a reporting parameter object
#' @export
imperfect_reporting_ind <- function(params) {
  
  # TODO - check for tau and eta?
  
  class(params) <- "imperfect_reporting_ind"
  return(params)
}

###############################################################################
#' Create a reporting graph from a social network
#'
#' Create a reporting graph with attributes added to vertices in the network
#' that have aggregate relational reports. Incorporate imperfect reporting
#' (i.e., false positives and negatives) at the individual level (varying
#' from person to person).
#' 
#' Note that this function relies upon the fact that the igraph object
#' \code{sim.graph} will have an attribute called \code{'sim.settings'},
#' which is a list with parameters describing the simulation.
#'
#' @details
#' For now, the only simulation parameter that makes a difference is
#' \eqn{\tau}. For this individual-level imperfect reporting, the
#' entry in the \code{reporting.params} list called \code{tau} should
#' be a function that takes two arguments: a vertex id and a graph.
#' It should return a value from 0 to 1, which is the tau for reports
#' from the given vertex.
#' 
#' @param reporting.params the reporting parameters
#' @param sim.graph the \code{igraph} object with the social network
#' @param stochastic if TRUE, then treat the reporting parameters as expected 
#' values from bernoulli trials for whether or not each edge is observed; otherwise,
#' reporting params are deterministic (but since edges are discrete, this could lead to
#' rounding issues -- ie, tau=0.8 for 3 edges would produce 2 edges observed)
#' @return the \code{igraph} object for the directed reporting graph 
#' @export
reporting_graph.imperfect_reporting_ind <- function(reporting.params, sim.graph, stochastic=FALSE) {
  
  
  # these should be functions that take the vertex and return values of
  # tau and eta, respectively
  tau <- reporting.params$tau
  eta <- reporting.params$eta
  
  rep.graph <- as.directed(sim.graph, mode='mutual')
  
  # figure out which vertices are in the hidden popn
  h.idx <- V(rep.graph)[in.H == 1] 
  f.idx <- V(rep.graph)[in.F == 1] 
  
  tolose <- c()
  
  # randomly remove edges 
  # NB: as.numeric(pot.edges) converts the igraph edgelist into a vector
  #     of edge ids, which we then sample
  for (f in f.idx) {
    
    # create a reporting graph by remove a fraction of the edges leading to
    # H from frame vertex f, according to the true positive rate tau
    pot.edges <- E(rep.graph)[f %->% h.idx]
    
    cur.tau <- tau(f, rep.graph)
    
    if(length(pot.edges) >= 1 & cur.tau < 1) {
      if (! stochastic) {
        # we'll take the fraction of edges as close as possible to the target (cur.tau)
        # but, because edges are discrete, this will often not be exact.
        # example: tau = 0.3 but there are two edges; here, we'll report one of the two edges,
        # so the realized tau is 0.5 and not 0.3
        tolose.idx <- sample(1:length(pot.edges), size=(1 - cur.tau)*length(pot.edges))
      } else {
        draws <- runif(length(pot.edges))
        tolose.idx <- which(draws > cur.tau)
      }
      tolose <- c(tolose, pot.edges[tolose.idx])
    }
  }
  
  # actually delete the edges
  rep.graph <- delete.edges(rep.graph, tolose)
  
  # TODO - implement eta
  
  class(rep.graph) <- c("igraph", class(sim.graph))
  
  # this returns a crazy vector that can't be right
  #tmp <- get.graph.attribute(rep.graph, 'groups')
  
  ## NB: BE CAREFUL HERE:
  ## counter-intuitively, for out-reports we want to use mode='in'
  ## and for in-reports we want to use mode='out'
  ## this cost a lot of time!
  ##
  # count up the reports in the reporting graph: 
  # ... out-reports (y)
  rep.graph <- report_edges(rep.graph, prefix='y.', mode="in")    
  
  # ... in-reports (v)
  rep.graph <- report_edges(rep.graph, prefix='v.', mode="out")    
  
  return(rep.graph)
  
}

