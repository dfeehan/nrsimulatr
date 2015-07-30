
###############################################################################
#' Create a set of parameters describing a census of F in graph
#'
#' This creates an object of class 'census_F'
#'
#' @param params a list of params, which must contain tau
#' @return a reporting parameter object
#' @export
frame_census <- function(params=NULL) {
    params <- list()
    class(params) <- "census_F"
    return(params)
}

###############################################################################
#' Get dataset from a census of F in reporting graph
#'
#' @param params sampling parameters; in this case, anything of the class 'census_F';
#'        its actual contents are ignored
#' @param reporting.graph the reporting graph
#' @return a dataset made from all of the nodes in F
#' @export
sample_graph.census_F <- function(params, reporting.graph) {

    rg.h <- induced.subgraph(reporting.graph, V(reporting.graph)[in.F == 1])

    rg.dat <- get.data.frame(rg.h, 'vertices')
    rg.dat$sampling.weight <- 1

    return(rg.dat)
}


###############################################################################
#' Create a set of parameters describing a census of H in graph
#'
#' This creates an object of class 'census_H'
#'
#' @param params a list of params, which must contain tau
#' @return a reporting parameter object
#' @export
hidden_census <- function(params=NULL) {
    params <- list()
    class(params) <- "census_H"
    return(params)
}

###############################################################################
#' Get dataset from a census of H in reporting graph
#'
#' @param params sampling parameters; in this case, anything of the class 'census_H';
#'        its actual contents are ignored
#' @param reporting.graph the reporting graph
#' @return a dataset made from all of the nodes in H
#' @export
sample_graph.census_H <- function(params, reporting.graph) {

    rg.h <- induced.subgraph(reporting.graph, V(reporting.graph)[in.H == 1])

    rg.dat <- get.data.frame(rg.h, 'vertices')
    rg.dat$sampling.weight <- 1

    return(rg.dat)
}

###############################################################################
#' Create a set of parameters describing a census of U in graph
#'
#' This creates an object of class 'census_U'
#'
#' @param params a list of params, which must contain tau
#' @return a reporting parameter object
#' @export
entire_census <- function(params=NULL) {
    params <- list()
    class(params) <- "census_U"
    return(params)
}

###############################################################################
#' Get dataset from a census of U in reporting graph
#'
#' @param params sampling parameters; in this case, anything of the class 'census_U';
#'        its actual contents are ignored
#' @param reporting.graph the reporting graph
#' @return a dataset made from all of the nodes in U
#' @export
sample_graph.census_U <- function(params, reporting.graph) {

    rg.dat <- get.data.frame(reporting.graph, 'vertices')
    rg.dat$sampling.weight <- 1

    return(rg.dat)
}


###############################################################################
#' Create a set of parameters describing a simple random sample of F in graph
#'
#' This creates an object of class 'srs_F', which has a simple
#' random sample without replacement of the frame population nodes
#' (NB: Sarndal et al call this SI sampling)
#'
#' @param params a list of params, which must contain sampling.frac, the sampling fraction
#' @return a reporting parameter object
#' @export
frame_srs <- function(params)
{

    # TODO - check that params has an entry sampling.frac

    class(params) <- "srs_F"
    return(params)
}

###############################################################################
#' Get dataset from a simple random sample of F in reporting graph
#'
#' @param params sampling parameters; see frame_srs
#' @param reporting.graph the reporting graph
#' @return a dataset made from a simple random sample of the nodes in F
#' @export
sample_graph.srs_F <- function(params, reporting.graph) {

    sampling.frac <- params$sampling.frac

    f.vertices <- V(reporting.graph)[in.F == 1]

    N.f <- length(f.vertices)
    sample.size <- floor(N.f*sampling.frac)

    sampled.idx <- sample(as.numeric(f.vertices), 
                          size=sample.size,
                          replace=FALSE)

    ## because the sample is a whole number, the empirical
    ## sampling fraction may be slightly different from the
    ## parameter
    empirical.frac <- sample.size / N.f

    rg.h <- induced.subgraph(reporting.graph, sampled.idx)

    rg.dat <- get.data.frame(rg.h, 'vertices')
    rg.dat$sampling.weight <- 1/empirical.frac

    return(rg.dat)
}

###############################################################################
#' Create a set of parameters describing a relative probability sample of H in graph
#'
#' This creates an object of class 'relprobsample_H', which has a relative
#' probability sample of H without replacement
#'
#' @param params a list of params, which must contain sampling.frac, the sampling fraction,
#'               and inclusion.trait, the trait to use for relative probabilities (typically degree)
#' @return a reporting parameter object
#' @export
hidden_relprob <- function(params) {

    # TODO - check that params has an entry sampling.frac

    class(params) <- "relprob_H"
    return(params)
}

###############################################################################
#' Get dataset from a relative probability sample of H in reporting graph
#'
#' @param params sampling parameters; see hidden_relprob
#' @param reporting.graph the reporting graph
#' @return a dataset made from a simple random sample of the nodes in F
#' @export
sample_graph.relprob_H <- function(params, reporting.graph) {

    sampling.frac <- params$sampling.frac
    inclusion.trait <- params$inclusion.trait

    h.vertices <- V(reporting.graph)[in.H == 1]

    trait.vals <- get.vertex.attribute(reporting.graph,
                                       inclusion.trait,
                                       h.vertices)

    trait.probs <- trait.vals / sum(trait.vals)

    N.h <- length(h.vertices)
    sample.size <- floor(N.h*sampling.frac)

    sampled.idx <- sample(as.numeric(h.vertices), 
                          size=sample.size,
                          replace=FALSE,
                          prob=trait.probs)

    ## because the sample is a whole number, the empirical
    ## sampling fraction may be slightly different from the
    ## parameter
    empirical.frac <- sample.size / N.h

    rg.h <- induced.subgraph(reporting.graph, sampled.idx)

    trait.sample <- get.vertex.attribute(rg.h,
                                         inclusion.trait)

    rg.dat <- get.data.frame(rg.h, 'vertices')
    rg.dat$sampling.weight <- 1/trait.sample

    return(rg.dat)
}




