
###############################################################################
#' draw a random graph from a stochastic blockmodel and label the vertices
#' 
#' Draw a random graph from a stochastic blockmodel using a k-vector
#' of block sizes and a k x k preference matrix
#'
#' The \code{igraph} object that gets returned is a graph drawn from
#' the stochastic block-model. The vertices each have a 'group' attribute
#' containing the name of the block they have been assigned to. These block
#' names default to 'group.1', 'group.2', ... if none are specified.
#'
#' @param block.sizes a vector with the number of vertices in each block
#' @param pref.matrix a matrix whose (i,j)th entry is the probability of
#'        an edge between a vertex in the block containing i and a vertex
#'        in the block containing j
#' @param block.names if not NULL, then names of the blocks
#' @return the \code{igraph} object with each vertex's group membership 
#'         given by a vertex attribute
#' @export
draw.sbm.graph <- function(block.sizes, pref.matrix, block.names=NULL) {

    # draw a random graph using a stochastic blockmodel
    g <- sbm.game(sum(block.sizes), pref.matrix, block.sizes)

    if (is.null(block.names)) {
        if (is.null(names(block.sizes))) {
           block.names <- paste0('group.', 1:length(block.sizes))
        } else {
            block.names <- names(block.sizes)
        }
    }

    # label vertices by group membership
    g <- label.vertices(g, block.sizes, block.names)

    return(g)

}

#####################################################################################
#' construct a 4-block preference matrix based on two independent categories
#' 
#' @description
#' TODO
#'
#' @param pi.within the probability of an edge between two people in the same group
#' @param rho the relative probability of an edge between two people who are not
#'        in the same group (see Description)
#' @param type (see Description)
#' @return a 4x4 preference matrix
pref.matrix.4group <- function(block.sizes, pi.within, rho, type,
                               inF=NULL, inH=NULL) {

    if (type=="simple") {
        
        # for the 'simple' type, we only care whether or not two blocks are exactly the same
        # the block preference matrix has pi.within*rho in all of its off-diagonal entries,
        # and pi.within on the diagonal
        pref.matrix <- matrix(pi.within*rho, 
                              nrow=length(block.sizes), ncol=length(block.sizes))
        diag(pref.matrix) <- pi.within

    } else if (type == "nested") {

        if (is.null(inF) | is.null(inH)) {
            stop("inF and inH must be speicifed to use the 'nested' option.")
        }

        # ties between F and not F are either from F to not F or from not F to F:
        Fmask <- (inF %*% t(1-inF)) + ((1-inF) %*% t(inF))
        Fmask <- 1 - Fmask * (1 - rho)

        # ties between H and not H are either from H to not H or from not H to H:
        Hmask <- (inH %*% t(1-inH)) + ((1-inH) %*% t(inH))
        Hmask <- 1 - Hmask * (1 - rho)

        # first make a matrix of the baseline prob of w/in group edge
        pref.matrix <- matrix(pi.within, nrow=length(block.sizes), ncol=length(block.sizes))
        # then apply the two masks which change probabilities of edges 
        # between F / not F and H / not H
        # (assuming these are independent)
        pref.matrix <- pref.matrix * Hmask * Fmask

    }

    colnames(pref.matrix) <- names(block.sizes)
    rownames(pref.matrix) <- names(block.sizes)

    return(pref.matrix)
}

###############################################################################
#' given parameter values, compute block sizes
#' 
#' given the parameter values N, p.H, and p.F, compute the number of vertices
#' in each of the four blocks: FnotH, FH, notFnotH, notFH
#' 
#'
#' @param sim.settings a list whose entries contain the parameter values
#'        which govern this random graph draw (see details)
#' @return a named vector; names are block names, values are block sizes
#' @export
block.sizes.4group <- function(sim.settings) {

    # there are more sophisticated ways of doing this, but
    # this has the advantage of being more readable, and
    # causing errors if parameters don't exist
    N <- sim.settings[['N']]
    p.F <- sim.settings[['p.F']]
    p.H <- sim.settings[['p.H']]
    pi.within <- sim.settings[['pi.within']]
    rho <- sim.settings[['rho']]

    # assume membership in F and H are independent
    N.FnotH <- floor(N * p.F * (1-p.H))
    N.notFnotH <- floor(N * (1-p.F) * (1-p.H))
    N.notFH <- floor(N * (1-p.F) * p.H)
    N.FH <- N - N.FnotH - N.notFnotH - N.notFH

    block.sizes <- c('FnotH' = N.FnotH,
                     'FH' = N.FH,
                     'notFnotH' = N.notFnotH,
                     'notFH' = N.notFH)

    return(block.sizes)

}


###############################################################################
#' draw a random graph from a stochastic blockmodel with four blocks
#' 
#' Draw a random graph from a stochastic blockmodel where the groups are
#' determined by membership in the frame population (F) and membership in
#' the hidden population (H). Membership in F is assumed to be
#' independent of membership in H, and vice-versa.
#'
#' @details
#' The four groups, as determined by membership in F and H, are
#' FH, FnotH, notF,notH, notFH.
#'
#' The \code{sim.settings} list should have entries
#' \itemize{
#' \item N the size of the entire population
#' \item p.F the prevalence of the frame population (N.F/N)
#' \item p.H the prevalence of the hidden population (N.H/N)
#' \item pi.within the probability of forming a tie for a pair of
#'       nodes in the same group
#' \item rho the relative probability of edge between two nodes that
#'       are not in the same group (for example, between a node in FnotH and FH).
#'       if rho is 0.8, then edges between FnotH and FH are 80% as likely as edges between
#'       FnotH and FnotH (or FH and FH)
#' \item tau the true positive rate
#' }
#'
#' The overall level of connectivity in the graph can be controlled by
#' varying pi.within. The extent to which members of the group are
#' homogenously mixed can be controlled by varying rho.
#' For example, when rho is 1, then all four groups are perfectly
#' mixed with the entire population. When rho is less than 1, members of the
#' same group tend to form ties with one another more than with others.
#' And when rho is greater than 1, members of the same group tend to
#' form ties with others more than with one another.
#'
#' @param sim.settings a list whose entries contain the parameter values
#'        which govern this random graph draw (see details)
#' @param type either 'simple' or 'nested', depending on the type of inhomogenous
#'        mixing; see \code{pref.matrix.4group}
#' @return an igraph graph, drawn according to the
#'         settings passed in
#' @export
draw.4group.graph <- function(sim.settings, type="simple") {

    # there are more sophisticated ways of doing this, but
    # this has the advantage of being more readable, and
    # causing errors if parameters don't exist
    pi.within <- sim.settings[['pi.within']]
    rho <- sim.settings[['rho']]

    block.sizes <- block.sizes.4group(sim.settings)

    pref.matrix <- pref.matrix.4group(block.sizes, 
                                      pi.within, rho, type,
                                      inF=c(1,1,0,0), inH=c(0,1,0,1))

    this.g <- draw.sbm.graph(block.sizes, pref.matrix)
    
    this.settings <- sim.settings
    this.res <- set.graph.attribute(this.g, "sim.settings", this.settings)

    this.res <- set.graph.attribute(this.res, "pref.matrix", pref.matrix)

    this.res <- set.graph.attribute(this.res, "block.sizes", block.sizes)

    V(this.res)$id <- 1:vcount(this.res)

    # report the graph edges
    this.res <- report.sbm.edges(this.res, 'd.')

    return(this.res)

}

###############################################################################
#' for each vertex, compute the number of reported connections to each block
#'
#' @details
#' 
#' This function is useful for computing edge counts (ie, reports) between
#' individual vertices and the blocks. These edge counts are the building
#' blocks of many network reporting estimators.
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

#####################################################################################
#' compute the expected number of connections (degree) between two groups under 4-block model
#' 
#' for groups A and B, this returns \eqn{d_{A,B} = d_{B,A}}. 
#'
#' @param pref.matrix the preference matrix
#' @param block.sizes vector with the number of vertices in each block
#' @param first.group vector with 0/1 entries indicating which blocks are in group A
#' @param second.group vector with 0/1 entries indicating which blocks are in group B
#' @return The expected number of edges between the two groups under the stochastic block 
#'         model given by pref.matrix. For example, 
expected.edgecount.4group <- function(pref.matrix, 
                                      block.sizes,
                                      first.group, 
                                      second.group) {

    fg <- first.group * block.sizes
    sg <- second.group * block.sizes
    both.groups <- (first.group * second.group * block.sizes)

    
    return(as.numeric(t(fg) %*% pref.matrix %*% sg -
                      t(both.groups)%*%diag(pref.matrix)))

}


#####################################################################################
#' compute the expected value for various quantities under four-group model
#' 
#' 
#'
#' @param these.params list containing the simulation parameters
#' @param type either "nested" or "simple"
#' @param inF vector with 1s for blocks in F, 0 otherwise
#' @param inH vector with 1s for blocks in H, 0 otherwise
#' @param inU vector with 1s for blocks in U, 0 otherwise
#' @return the expected value of various quantities
#' @export
ev.4group <- function(these.params, 
                      type="nested",
                      inF=c(1,1,0,0),
                      inH=c(0,1,0,1),
                      inU=c(1,1,1,1)) {

    these.block.sizes <- block.sizes.4group(these.params)

    this.pref.matrix <- pref.matrix.4group(these.block.sizes,
                                           these.params$pi.within,
                                           these.params$rho,
                                           type=type,
                                           inF=this.inF,
                                           inH=this.inH)

    res <- these.params

    res$dbar.F.F <- ev.4group.mean(this.pref.matrix, these.block.sizes,
                                   inF, inF)
    res$dbar.notF.F <- ev.4group.mean(this.pref.matrix, these.block.sizes,
                                   1-inF, inF)
    res$d.F.H <- ev.4group.tot(this.pref.matrix, these.block.sizes,
                               inF, inH)
    res$dbar.U.F <- ev.4group.mean(this.pref.matrix, these.block.sizes,
                                   inU, inF)
    res$dbar.F.U <- ev.4group.mean(this.pref.matrix, these.block.sizes,
                                   inF, inU)
    res$dbar.U.U <- ev.4group.mean(this.pref.matrix, these.block.sizes,
                                   inU, inU)
    res$dbar.H.U <- ev.4group.mean(this.pref.matrix, these.block.sizes,
                                   inH, inU)
    res$dbar.H.F <- ev.4group.mean(this.pref.matrix, these.block.sizes,
                                   inH, inF)

    res$phi.F <- res$dbar.F.F / res$dbar.U.F
    res$delta.F <- res$dbar.H.F / res$dbar.F.F
    res$deltaXphi.F <- res$phi.F * res$delta.F
    res$total <- res$phi.F * res$delta.F * res$tau

    return(data.frame(res))

}

#####################################################################################
#' compute the expected value for various quantities under four-group model
#' 
#' for groups A and B, these functions return \eqn{d_{A,B}} and \eqn{\bar{d}_{A,B}}
#'
#' @param pm the preference matrix
#' @param bs the block sizes
#' @param fg 0/1 indicator vector for blocks in the first index (A in the description above)
#' @param sg 0/1 indicator vector for blocks in the second index (B in the description above)
#' @return the expected value of the mean or the total number of connections
#'         from the first group to the second. Means are taken with respect to the members
#'         of the first index
ev.4group.tot <- function(pm, bs, fg, sg) {
    return(expected.edgecount.4group(pm, bs, fg, sg))
}

#' @rdname ev.4group.tot
ev.4group.mean <- function(pm, bs, fg, sg) {
    return(expected.edgecount.4group(pm, bs, fg, sg)/sum(bs*fg))
}

