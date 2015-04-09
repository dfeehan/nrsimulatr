
## TODO - comment this
## TODO - move the checks to the -test file

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
#' @param sim.graph the \code{igraph} object with the social network
#' @param hidden.popn the name(s) of blocks containing hidden population members
#' @param frame.popn the name(s) of blocks containing frame population members
#' @return the \code{igraph} object with each vertex's group membership 
#'         given by vertex label
#' @export
reporting.graph <- function(sim.graph, hidden.popn, frame.popn) {

    sim.settings <- get.graph.attribute(sim.graph, 'sim.settings')
    tau <- sim.settings[['tau']]

    rep.graph <- as.directed(sim.graph, mode='mutual')

    # figure out which vertices are in the hidden popn
    h.idx <- V(rep.graph)[group %in% hidden.popn] 
    f.idx <- V(rep.graph)[group %in% frame.popn] 

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

    ## NB: BE CAREFUL HERE:
    ## counter-intuitively, for out-reports we want to use mode='in'
    ## and for in-reports we want to use mode='out'
    ## this cost a lot of time!
    ##
    # count up the reports in the reporting graph: 
    # ... out-reports (y)
    rep.graph <- report.sbm.edges(rep.graph, prefix='y.', mode="in")    
    # ... in-reports (v)
    rep.graph <- report.sbm.edges(rep.graph, prefix='v.', mode="out")    

    return(rep.graph)

}

