###############################################################################
#' compute network reporting estimates from a randomly simulated graph
#'
#' note that this code relies upon the fact that the vertices should have
#' two 0/1 attributes, in.H and in.F, which indicate whether or not each one
#' is in the hidden population and the frame population
#'
#' @param rg.dat dataset from a reporting graph object
#' @param weights weights that account for the sampling (if any);
#'        this should be a vector that is the same length as the number of rows in rg.dat
#' @return a data frame with estimates and several related quantities
#' @export
sbm_estimates <- function(data, weights) {

    rg.dat <- get.data.frame(reporting.graph, 'vertices')

    true.N.H <- sum(rg.dat$in.H)

    nsum.ests <- rg.dat %>% 
                 summarise(tot.y = sum(in.F*(y.FH + y.notFH)),
                           dbar.F.F = sum(in.F*(d.FH + d.FnotH))/sum(in.F),
                           dbar.U.F = mean(d.FH + d.FnotH),
                           dbar.H.F = sum(in.H*(d.FH + d.FnotH))/sum(in.H),
                           vbar.H.F = sum(in.F*(y.FH + y.notFH))/sum(in.H),
                           dbar.FH = sum(in.F*in.H*d.degree)/sum(in.F*in.H),
                           dbar.FnotH = sum(in.F*(1-in.H)*d.degree)/sum(in.F*(1-in.H)),
                           dbar.notFnotH = sum((1-in.F)*(1-in.H)*d.degree)/
                                           sum((1-in.F)*(1-in.H)),
                           dbar.notFH = sum((1-in.F)*in.H*d.degree)/sum((1-in.F)*in.H),
                           dbar = mean(d.degree)) %>%
                ## compute the adjustment factors, as well as the basic
                ## and generalized estimates
                mutate(phi.F = dbar.F.F / dbar.U.F,
                       delta.F = dbar.H.F / dbar.F.F,
                       tau.F = vbar.H.F / dbar.H.F,
                       basic.est = tot.y / dbar.U.F,
                       generalized.est = tot.y / vbar.H.F,
                       ## compute the bias in the basic estimator
                       basic.bias = basic.est - true.N.H,
                       ## and double check that our identity works
                       id.check = basic.est * (1/phi.F) * (1/delta.F) * (1/tau.F))

    return(nsum.ests)

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
report_edges.sbm <- function(g, prefix='d.', mode="all") {

    # compute the degree of each vertex
    # note that for undirected graphs, mode='out' has no effect; this just
    # returns the degree
    #V(g)$degree <- degree(g, mode="out")
    g <- set.vertex.attribute(g, paste0(prefix, 'degree'), value=degree(g, mode=mode))

    # and compute a matrix of indicator variables with the group membership of each vertex
    gp.lookup <- get.data.frame(g, 'vertices')

    #block.names <- unique(paste(gp.lookup$group))
    block.names <- names(graph.attributes(g)$block.sizes)

    gp.lookup$group <- factor(gp.lookup$group, levels=block.names)

    gp.lookup$id <- 1:nrow(gp.lookup)

    # (make this wide, so that we can just sum columns to get totals)
    gp.lookup <- gp.lookup %>% 
                 select(id, group) %>%
                 mutate(present=1) %>% 
                 spread(group, present, fill=0, drop=FALSE) %>% 
                 arrange(id)

    gp.lookup <- as.matrix(gp.lookup[,block.names])

    # go through each vertex, and add 1 to its neighbors' counts of
    # connections based on which group it's in

    ## get reports (out-edges)
    # start with a matrix of 0s for each vertex's connection to each other group
    ubertally <- matrix(0, nrow=length(V(g)), ncol=length(block.names))
    colnames(ubertally) <- block.names

    for(i in 1:length(V(g))) {
        these.nei <- V(g)[nei(i, mode=mode)]
        ubertally[these.nei,] <- t(t(ubertally[these.nei,]) + gp.lookup[i,])
    }

    colnames(ubertally) <- paste0(prefix, colnames(ubertally))

    # add each vertex's number of connections as a vertex attribute called,
    # eg, "d.group1", "d.group2", etc...
    for(this.name in colnames(ubertally)) {
        g <- set.vertex.attribute(g, this.name, value=ubertally[,this.name])
    }

    return(g)

}

################## TEMP - BELOW HERE IS OLD STUFF...




###############################################################################
#' for each vertex, compute the number of reported connections to each block
#'
#' @details
#' 
#' This function is useful for computing edge counts (ie, reports) between
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
report.sbm.edges <- function(g, prefix='d.', mode="all") {

    # compute the degree of each vertex
    # note that for undirected graphs, mode='out' has no effect; this just
    # returns the degree
    #V(g)$degree <- degree(g, mode="out")
    g <- set.vertex.attribute(g, paste0(prefix, 'degree'), value=degree(g, mode=mode))

    # and compute a matrix of indicator variables with the group membership of each vertex
    gp.lookup <- get.data.frame(g, 'vertices')

    #block.names <- unique(paste(gp.lookup$group))
    block.names <- names(graph.attributes(g)$block.sizes)

    gp.lookup$group <- factor(gp.lookup$group, levels=block.names)

    gp.lookup$id <- 1:nrow(gp.lookup)

    # (make this wide, so that we can just sum columns to get totals)
    gp.lookup <- gp.lookup %>% 
                 select(id, group) %>%
                 mutate(present=1) %>% 
                 spread(group, present, fill=0, drop=FALSE) %>% 
                 arrange(id)

    gp.lookup <- as.matrix(gp.lookup[,block.names])

    # go through each vertex, and add 1 to its neighbors' counts of
    # connections based on which group it's in

    ## get reports (out-edges)
    # start with a matrix of 0s for each vertex's connection to each other group
    ubertally <- matrix(0, nrow=length(V(g)), ncol=length(block.names))
    colnames(ubertally) <- block.names

    for(i in 1:length(V(g))) {
        these.nei <- V(g)[nei(i, mode=mode)]
        ubertally[these.nei,] <- t(t(ubertally[these.nei,]) + gp.lookup[i,])
    }

    colnames(ubertally) <- paste0(prefix, colnames(ubertally))

    for(this.name in colnames(ubertally)) {
        g <- set.vertex.attribute(g, this.name, value=ubertally[,this.name])
    }

    return(g)

}

###############################################################################
#' compute network reporting estimates from a randomly simulated graph
#'
#' @param reporting.graph the reporting graph object
#' @param hidden.popn the names of block(s) containing the hidden population
#' @param frame.popn the names of block(s) containing the frame population
#' @return a data frame with estimates and several related quantities
#' @export
reporting.estimates <- function(reporting.graph, hidden.popn, frame.popn) {

    rg.dat <- get.data.frame(reporting.graph, 'vertices')
    rg.dat <- rg.dat %>% mutate(in.H = as.numeric(group %in% hidden.popn),
                                in.F = as.numeric(group %in% frame.popn))

    true.N.H <- sum(rg.dat$in.H)

    nsum.ests <- rg.dat %>% 
                 summarise(tot.y = sum(in.F*(y.FH + y.notFH)),
                           dbar.F.F = sum(in.F*(d.FH + d.FnotH))/sum(in.F),
                           dbar.U.F = mean(d.FH + d.FnotH),
                           dbar.H.F = sum(in.H*(d.FH + d.FnotH))/sum(in.H),
                           vbar.H.F = sum(in.F*(y.FH + y.notFH))/sum(in.H),
                           dbar.FH = sum(in.F*in.H*d.degree)/sum(in.F*in.H),
                           dbar.FnotH = sum(in.F*(1-in.H)*d.degree)/sum(in.F*(1-in.H)),
                           dbar.notFnotH = sum((1-in.F)*(1-in.H)*d.degree)/
                                           sum((1-in.F)*(1-in.H)),
                           dbar.notFH = sum((1-in.F)*in.H*d.degree)/sum((1-in.F)*in.H),
                           dbar = mean(d.degree)) %>%
                ## compute the adjustment factors, as well as the basic
                ## and generalized estimates
                mutate(phi.F = dbar.F.F / dbar.U.F,
                       delta.F = dbar.H.F / dbar.F.F,
                       tau.F = vbar.H.F / dbar.H.F,
                       basic.est = tot.y / dbar.U.F,
                       generalized.est = tot.y / vbar.H.F,
                       ## compute the bias in the basic estimator
                       basic.bias = basic.est - true.N.H,
                       ## and double check that our identity works
                       id.check = basic.est * (1/phi.F) * (1/delta.F) * (1/tau.F))

    return(nsum.ests)

}

