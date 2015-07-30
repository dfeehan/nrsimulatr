
###############################################################################
#' Create a set of parameters describing reporting behavior that is not perfectly accurate
#'
#' This creates an object of class 'imperfect_reporting' whose parameter values
#' correspond to perfect reporting.
#'
#' @details
#' For now, the only simulation parameter that makes a difference is
#' \eqn{\tau}.
#'
#' @param params a list of params, which must contain tau
#' @return a reporting parameter object
#' @export
perfect_reporting <- function(params=NULL) {
    params <- list(tau=1, eta=1)
    class(params) <- "imperfect_reporting"
    return(params)
}

###############################################################################
#' Create a set of parameters describing reporting behavior that is not perfectly accurate
#'
#'
#' @details
#' For now, the only simulation parameter that makes a difference is
#' \eqn{\tau}.
#'
#' @param params a list of params, which must contain tau
#' @return a reporting parameter object
#' @export
imperfect_reporting <- function(params) {

    # TODO - check for tau and eta?

    class(params) <- "imperfect_reporting"
    return(params)
}

###############################################################################
#' Create a reporting graph from a social network
#'
#' Note that this function relies upon the fact that the igraph object
#' \code{sim.graph} will have an attribute called \code{'sim.settings'},
#' which is a list with parameters describing the simulation.
#'
#' @details
#' For now, the only simulation parameter that makes a difference is
#' \eqn{\tau}.
#'
#' @param reporting.params the reporting parameters
#' @param sim.graph the \code{igraph} object with the social network
#' @return the \code{igraph} object for the directed reporting graph 
#' @export
reporting_graph.imperfect_reporting <- function(reporting.params, sim.graph) {

    tau <- reporting.params$tau
    eta <- reporting.params$eta

    rep.graph <- as.directed(sim.graph, mode='mutual')


    # figure out which vertices are in the hidden popn
    h.idx <- V(rep.graph)[in.H == 1] 
    f.idx <- V(rep.graph)[in.F == 1] 

    # randomly remove edges 
    # NB: as.numeric(pot.edges) converts the igraph edgelist into a vector
    #     of edge ids, which we then sample
    if (tau < 1) {
      # create a reporting graph by remove a fraction of the edges leading to
      # H from F, according to the true positive rate tau
      pot.edges <- E(rep.graph)[f.idx %->% h.idx]

      tolose <- sample(as.numeric(pot.edges), size=(1 - tau)*length(pot.edges))
      rep.graph <- delete.edges(rep.graph, tolose)
    }

    # TODO - implement eta

    class(rep.graph) <- c("igraph", class(sim.graph))

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

    ## TODO - decide what to do about what class the reporting
    ## graph is

    return(rep.graph)

}






