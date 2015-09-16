###############################################################################
#' create a set of simulation parameters for the stochastic block model
#' 
#' @param params the simulation parameters
#' @param type currently, this can be one of '4group_1param_simple', '4group_1param_nested', or '4group_2param'; see the vignette for more info
#' @return an object with the simulation parameters
#' @export
sbm_params <- function(params, type) {
    
    # TODO - eventually, check that the params are all there?

    class(params) <- c(type, "sbm_4group", "sbm")
    block.sizes <- sbm_block_sizes(params)
    params <- c(params, block.sizes)

    class(params) <- c(type, "sbm_4group", "sbm")
    pref.matrix <- sbm_pref_matrix(params)
    params$pref.matrix <- pref.matrix

    return(params)
}

###############################################################################
#' given parameter values for a stochastic block model, compute block sizes (generic)
#' 
#' given the parameter values of a stochastic block model, return a named
#' list of the block sizes
#'
#' @param params a list whose entries contain the parameter values
#'        which govern this random graph draw (see details)
#' @return a named vector; names are block names, values are block sizes
#' @export
sbm_block_sizes <- function(params) {
    UseMethod("sbm_block_sizes")
}

###############################################################################
#' given parameter values, compute block sizes
#' 
#' given the parameter values N, p.H, and p.F, compute the number of vertices
#' in each of the four blocks: FnotH, FH, notFnotH, notFH
#'
#' @param params a list whose entries contain the parameter values
#'        which govern this random graph draw (see details)
#' @return a list; block.sizes has a vector whose names are block names, values are block sizes;
#'         gps.in.F is a vector with the names of groups whose members are in F; and
#'         gps.in.H is a vector with the names of groups whose members are in H
#' @export
sbm_block_sizes.sbm_4group <- function(params) {

    # there are more sophisticated ways of doing this, but
    # this has the advantage of being more readable, and
    # causing errors if parameters don't exist
    N <- params[['N']]
    p.F <- params[['p.F']]
    p.H <- params[['p.H']]

    # if the parameters 'p.F.given.H' is either not there (NULL) or NA
    # then assume independent membership in the two groups
    p.F.given.H <- ifelse(! (is.null(params[['p.F.given.H']])||is.na(params[['p.F.given.H']])),
                          params[['p.F.given.H']],
                          p.F)

    zeta <- params[['zeta']]
    rho <- params[['rho']]

    # assume membership in F and H are independent
    N.F <- floor(N*p.F)
    N.H <- floor(N*p.H)

    # if p.F.given.H * N.H > N.F, we only take
    # N.F ppl
    N.FH <- floor(min(p.F.given.H*N.H, N.F))
    N.FnotH <- N.F - N.FH
    N.notFH <- N.H - N.FH
    N.notFnotH <- N - N.FH - N.FnotH - N.notFH

    #N.FnotH <- floor(N * p.F * (1-p.H))
    #N.notFnotH <- floor(N * (1-p.F) * (1-p.H))
    #N.notFH <- floor(N * (1-p.F) * p.H)
    #N.FH <- N - N.FnotH - N.notFnotH - N.notFH

    block.sizes <- c('FnotH' = N.FnotH,
                     'FH' = N.FH,
                     'notFnotH' = N.notFnotH,
                     'notFH' = N.notFH)

    # get the name of the groups in F
    gps.in.F <- c('FnotH', 'FH')

    # get the name of the groups in H
    gps.in.H <- c('FH', 'notFH')

    return(list(block.sizes=block.sizes,
                gps.in.F=gps.in.F,
                gps.in.H=gps.in.H))

}

#####################################################################################
#' construct a preference matrix for a stochastic block model
#' 
#' @description
#' TODO
#'
#' @param params the parameters specific to this type of model
#' @return a preference matrix
sbm_pref_matrix <- function(params) {
    UseMethod("sbm_pref_matrix")
}

#####################################################################################
#' construct a 4-block preference matrix based on two independent categories, simple
#' 
#' @description
#' TODO
#'
#' @param list of parameters, including block.sizes, zeta and rho
#' @return a 4x4 preference matrix
sbm_pref_matrix.4group_1param_simple <- function(params) {

    zeta <- params$zeta
    rho <- params$rho
    block.sizes <- params$block.sizes

    pref.matrix <- matrix(zeta*rho, 
                          nrow=length(block.sizes), ncol=length(block.sizes))
    diag(pref.matrix) <- zeta

    colnames(pref.matrix) <- names(block.sizes)
    rownames(pref.matrix) <- names(block.sizes)

    return(pref.matrix)
}

#####################################################################################
#' construct a 4-block preference matrix based on two independent categories, nested
#' 
#' @description
#' TODO
#'
#' @param list of parameters, including block.sizes, zeta, rho, gps.in.F, and gps.in.H
#' @return a 4x4 preference matrix
sbm_pref_matrix.4group_1param_nested <- function(params) {

    zeta <- params$zeta
    block.sizes <- params$block.sizes
    rho <- params$rho

    gp.in.F <- as.numeric(names(block.sizes) %in% params$gps.in.F)
    gp.in.H <- as.numeric(names(block.sizes) %in% params$gps.in.H)

    # ties between F and not F are either from F to not F or from not F to F:
    Fmask <- (gp.in.F %*% t(1-gp.in.F)) + ((1-gp.in.F) %*% t(gp.in.F))
    Fmask <- 1 - Fmask * (1 - rho)

    # ties between H and not H are either from H to not H or from not H to H:
    Hmask <- (gp.in.H %*% t(1-gp.in.H)) + ((1-gp.in.H) %*% t(gp.in.H))
    Hmask <- 1 - Hmask * (1 - rho)

    # first make a matrix of the baseline prob of w/in group edge
    pref.matrix <- matrix(zeta, nrow=length(block.sizes), ncol=length(block.sizes))
    # then apply the two masks which change probabilities of edges 
    # between F / not F and H / not H
    # (assuming these are independent)
    pref.matrix <- pref.matrix * Hmask * Fmask

    colnames(pref.matrix) <- names(block.sizes)
    rownames(pref.matrix) <- names(block.sizes)

    return(pref.matrix)
}

#####################################################################################
#' construct a 4-block preference matrix based on two independent categories,
#' and two independent interaction parameters
#' 
#' @description
#' TODO
#'
#' @param list of parameters, including block.sizes, zeta, rho, xi, gps.in.F, and gps.in.H
#' @return a 4x4 preference matrix
sbm_pref_matrix.4group_2param <- function(params) {

    zeta <- params$zeta
    block.sizes <- params$block.sizes
    # rho is the penalty for prob of edge between F vs not F
    rho <- params$rho
    # xi is the penalty for prob of edge between H vs not H
    xi <- params$xi

    gp.in.F <- as.numeric(names(block.sizes) %in% params$gps.in.F)
    gp.in.H <- as.numeric(names(block.sizes) %in% params$gps.in.H)

    # ties between F and not F are either from F to not F or from not F to F:
    Fmask <- (gp.in.F %*% t(1-gp.in.F)) + ((1-gp.in.F) %*% t(gp.in.F))
    Fmask <- 1 - Fmask * (1 - xi)

    # ties between H and not H are either from H to not H or from not H to H:
    Hmask <- (gp.in.H %*% t(1-gp.in.H)) + ((1-gp.in.H) %*% t(gp.in.H))
    Hmask <- 1 - Hmask * (1 - rho)

    # first make a matrix of the baseline prob of w/in group edge
    pref.matrix <- matrix(zeta, nrow=length(block.sizes), ncol=length(block.sizes))
    # then apply the two masks which change probabilities of edges 
    # between F / not F and H / not H
    # (assuming these are independent)
    pref.matrix <- pref.matrix * Hmask * Fmask

    colnames(pref.matrix) <- names(block.sizes)
    rownames(pref.matrix) <- names(block.sizes)

    return(pref.matrix)
}

###############################################################################
#' draw a random graph from a 4group_oneparam_simple stochastic blockmodel
#' 
#' Draw a random graph from a stochastic blockmodel using a k-vector
#' of block sizes and a k x k preference matrix
#'
#' This routine is agnostic about the algorithm that was used to produce
#' the preference matrix and the block sizes, so it could potentially be
#' used with a lot of different models.
#'
#' The \code{igraph} object that gets returned is a graph drawn from
#' the stochastic block-model. The vertices each have a 'group' attribute
#' containing the name of the block they have been assigned to. These block
#' names default to 'group.1', 'group.2', ... if none are specified.
#'
#' @param params a list containing block.sizes, pref.matrix, and optionally block.names
#' @return the \code{igraph} object with each vertex's group membership 
#'         given by a vertex attribute
#' @export
generate_graph.sbm <- function(params) {

    block.sizes <- params$block.sizes
    pref.matrix <- params$pref.matrix
    block.names <- params$block.names
    gps.in.F <- params$gps.in.F
    gps.in.H <- params$gps.in.H

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
    g <- sbm_label_vertices(g, block.sizes, block.names)

    # label vertices w/ dummy variables for membership in F
    V(g)$in.F <- 0
    V(g)[V(g)$group %in% gps.in.F]$in.F <- 1

    # ... and for membership in H
    V(g)$in.H <- 0
    V(g)[V(g)$group %in% gps.in.H]$in.H <- 1

    class(g) <- append(class(g), class(params))

    return(g)

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
#' \item p.F.given.H the probability of being in F, given that a node is in H;
#'       if not specified, assume the probabilities are independent, so that
#'       p.F.given.H = p.F
#' \item zeta the probability of forming a tie for a pair of
#'       nodes in the same group
#' \item rho the relative probability of edge between two nodes that
#'       are not in the same group (for example, between a node in FnotH and FH).
#'       if rho is 0.8, then edges between FnotH and FH are 80% as likely as edges between
#'       FnotH and FnotH (or FH and FH)
#' \item tau the true positive rate
#' }
#'
#' The overall level of connectivity in the graph can be controlled by
#' varying zeta. The extent to which members of the group are
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
#'        mixing; see \code{pref.matrix.4group.1param}
#' @return an igraph graph, drawn according to the
#'         settings passed in
#' @export
generate_graph.sbm_4group <- function(params) {

    this.g <- generate_graph.sbm(params)

    this.res <- set.graph.attribute(this.g, "params", params)
    this.res <- set.graph.attribute(this.res, "pref.matrix", params$pref.matrix)
    this.res <- set.graph.attribute(this.res, "block.sizes", params$block.sizes)

    V(this.res)$id <- 1:vcount(this.res)

    # report the graph edges
    this.res <- report_edges(this.res, 'd.')

    return(this.res)

}

#' assign labels to vertices in a graph
#'
#' @param graph the \code{igraph} graph
#' @param block.sizes a vector with the number of nodes in each group
#' @param block.names a vector (same length as \code{block.sizes}) with names of each group
#' @export
sbm_label_vertices <- function(graph, block.sizes, block.names) {

    V(graph)$group <- "none"

    break.start <- c(0, cumsum(block.sizes))
    break.ends <- cumsum(block.sizes)

    break.start <- break.start[-length(break.start)] + 1

    for (i in 1:length(block.sizes)) {
        # when this condition is FALSE, block i has 0 members (so we don't label anything)
        if(break.start[i] <= break.ends[i]) {
            V(graph)[break.start[i]:break.ends[i]]$group <- block.names[i]
        }
    }

    return(graph)

}

#####################################################################################
#' generic function for getting expected edgecounts
#' 
#' @param params parameter vector describing the model (the class of this param vector
#'        determines which method gets run)
expected_edgecount <- function(params, ...) {
    UseMethod("expected_edgecount")
}

#####################################################################################
#' compute the expected number of connections (degree) between two groups under 4-block model
#' 
#' for groups A and B, this returns \eqn{d_{A,B} = d_{B,A}}. 
#'
#' it does not matter how the preference matrix and block sizes were generated,
#' so this function could be used with many different models.
#'
#' @param params a parameter vector, which must have pref.matrix and block.sizes
#' @param first.group vector with the names of the block that blocks that are in group A
#' @param second.group vector with the names of the blocks that are in group B
#' @return The expected number of edges between the two groups under the stochastic block 
#'         model given by pref.matrix. For example, 
expected_edgecount.sbm_4group <- function(params, 
                                          first.group, 
                                          second.group) {

    pref.matrix <- params$pref.matrix
    block.sizes <- params$block.sizes

    fg.idx <- as.numeric(names(params$block.sizes) %in% first.group)
    sg.idx <- as.numeric(names(params$block.sizes) %in% second.group)

    fg <- fg.idx * block.sizes
    sg <- sg.idx * block.sizes
    both.groups <- (fg.idx * sg.idx * block.sizes)

    return(as.numeric(t(fg) %*% pref.matrix %*% sg -
                      t((both.groups^2 + both.groups)/2)%*%diag(pref.matrix)))

}

#####################################################################################
#' generic fn for getting expected values from stochastic block model
#' 
#' @param params the parameter object
sbm_ev <- function(params, ...) {
    UseMethod("sbm_ev")
}

#####################################################################################
#' compute the expected value for various quantities under four-group model
#' 
#' 
#'
#' @param params list/object containing the simulation parameters
#' @param inF vector with 1s for blocks in F, 0 otherwise
#' @param inH vector with 1s for blocks in H, 0 otherwise
#' @param inU vector with 1s for blocks in U, 0 otherwise
#' @return the expected value of various quantities
#' @export
sbm_ev.sbm_4group <- function(params, 
                              inF=c('FnotH', 'FH'),
                              inH=c('FH', 'notFH'),
                              inU=c('FnotH', 'FH', 'notFnotH', 'notFH')) {


    these.block.sizes <- params$block.sizes
    this.pref.matrix <- params$pref.matrix

    not.inF <- inU[ ! inF %in% inU ]

    N.F <- sum(these.block.sizes[inF])
    N.H <- sum(these.block.sizes[inH])
    N <- sum(these.block.sizes)

    res <- params[c('N', 'p.F', 'p.H', 'p.F.given.H', 'zeta', 'rho', 'tau')]

    res$xi <- ifelse(! is.null(params$xi), params$xi, NA)

    res$d.F.H <- expected_edgecount(params, inF, inH)

    # averages are always wrt the first subscript
    res$dbar.F.F <- expected_edgecount(params, inF, inF) / N.F
    res$dbar.notF.F <- expected_edgecount(params, not.inF, inF) / (N - N.F)
    res$dbar.U.F <- expected_edgecount(params, inU, inF) / N
    res$dbar.F.U <- expected_edgecount(params, inF, inU) / N.F
    res$dbar.U.U <- expected_edgecount(params, inU, inU) / N
    res$dbar.H.U <- expected_edgecount(params, inH, inU) / N.H
    res$dbar.H.F <- expected_edgecount(params, inH, inF) / N.H

    res$phi.F <- res$dbar.F.F / res$dbar.U.F
    res$delta.F <- res$dbar.H.F / res$dbar.F.F
    res$deltaXphi.F <- res$phi.F * res$delta.F
    res$total <- res$phi.F * res$delta.F * res$tau

    return(data.frame(res))

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

    block.names <- names(graph.attributes(g)$block.sizes)

    adj.matrix <- get.adjacency(g)

    vertex.groups <- get.vertex.attribute(g, 'group')

    ## go through and get reported connections to each block
    for (b in block.names) {

        in.b <- as.numeric(vertex.groups == b)

        res.name <- paste0(prefix, b)

        if (mode=='out') {
          # in this case, we want visibilities
          # which are A^T x; we'll compute them via
          # we'll compute x^T A instead
          num.rep.conns <- as.numeric(t(in.b) %*% adj.matrix)
        } else {
          num.rep.conns <- as.numeric(adj.matrix %*% in.b)
        }

        g <- set.vertex.attribute(g, res.name,
                                  value=num.rep.conns)
    }

    return(g)

}

###############################################################################
#' compute network reporting estimates from a randomly simulated stochastic block model
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
sbm_census_quantities <- function(data, weights=data[,'sampling.weight']) {

    rg.dat <- data
    rg.dat$.weights <- weights

    true.N.H <- sum(rg.dat$in.H * rg.dat$.weights)

    nsum.ests <- rg.dat %>% 
                 summarise(tot.y = sum(.weights*in.F*(y.FH + y.notFH)),
                           dbar.F.F = sum(.weights*in.F*(d.FH + d.FnotH))/sum(.weights*in.F),
                           dbar.U.F = sum(.weights*(d.FH + d.FnotH))/sum(.weights),
                           dbar.H.F = sum(.weights*in.H*(d.FH + d.FnotH))/sum(.weights*in.H),
                           vbar.H.F = sum(.weights*in.F*(y.FH + y.notFH))/sum(.weights*in.H),
                           dbar.FH = sum(.weights*in.F*in.H*d.degree)/sum(.weights*in.F*in.H),
                           dbar.FnotH = sum(.weights*in.F*(1-in.H)*d.degree)/
                                        sum(.weights*in.F*(1-in.H)),
                           dbar.notFnotH = sum(.weights*(1-in.F)*(1-in.H)*d.degree)/
                                           sum(.weights*(1-in.F)*(1-in.H)),
                           dbar.notFH = sum(.weights*(1-in.F)*in.H*d.degree)/
                                        sum(.weights*(1-in.F)*in.H),
                           dbar = sum(.weights*d.degree)/sum(.weights)) %>%
                ## compute the adjustment factors, as well as the basic
                ## and generalized estimates
                mutate(phi.F = dbar.F.F / dbar.U.F,
                       delta.F = dbar.H.F / dbar.F.F,
                       tau.F = vbar.H.F / dbar.H.F,
                       census.basic.est = tot.y / dbar.U.F,
                       census.adapted.est = tot.y / dbar.F.F,
                       census.generalized.est = tot.y / vbar.H.F,
                       ## compute the bias in the basic estimator
                       census.basic.bias = census.basic.est - true.N.H,
                       ## and double check that our identity works
                       id.check = census.basic.est * (1/phi.F) * (1/delta.F) * (1/tau.F))

    return(nsum.ests)

}


