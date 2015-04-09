
OLD.draw.sbm.graph <- function(block.sizes, pref.matrix, block.names=NULL) {

    # draw a random graph using a stochastic blockmodel
    g <- sbm.game(N, pref.matrix, block.sizes)

    if (is.null(block.names)) {
        if (is.null(names(block.sizes))) {
           block.names <- paste0('group.', 1:length(block.sizes))
        } else {
            block.names <- names(block.sizes)
        }
    }

    # label vertices by group membership
    g <- label.vertices(g, block.sizes, block.names)

    # compute the degree of each vertex
    V(g)$degree <- degree(g)

    # and compute a matrix of indicator variables with the group membership of each vertex
    gp.lookup <- get.data.frame(g, 'vertices')
    gp.lookup$id <- 1:nrow(gp.lookup)
    # (make this wide, so that we can just sum columns to get totals)
    gp.lookup <- gp.lookup %>% mutate(present=1) %>% spread(group, present, fill=0) %>% arrange(id)

    gp.lookup <- as.matrix(gp.lookup[,block.names])

    ## this was the first draft -- it appears to be much slower than the version below
    #for(i in 1:length(V(g))) {
    #
    #    # get i's neighbors
    #    if (V(g)[i]$degree > 0) {
    #        these.nbrs <- V(g)[nei(i)]
    #        tally <- gp.lookup %>% filter(id %in% these.nbrs)
    #        tally <- colSums(tally %>% select(-degree, -id))
    #    } else {
    #        # we have an isolate
    #        tally <- rep(0, length(block.names))
    #        names(tally) <- block.names
    #    }
    #
    #    # add the number of connections to each vertex's attributes
    #    for(this.name in block.names) {
    #        g <- set.vertex.attribute(g, name=paste0('y.', this.name), index=i, value=tally[this.name])
    #    }
    #
    #}

    ## second version: go through each vertex, and add 1 to its neighbors counts of
    ## connections based on which group it's in

    # start with a matrix of 0s for each vertex's connection to each other group
    ubertally <- matrix(0, nrow=length(V(g)), ncol=length(block.names))
    colnames(ubertally) <- block.names

    #for(i in 1:length(V(g))) {
    for(i in 1:length(V(g))) {
        these.nei <- V(g)[nei(i)]
        ubertally[these.nei,] <- t(t(ubertally[these.nei,]) + gp.lookup[i,])
    }

    colnames(ubertally) <- paste0('d.', colnames(ubertally))

    for(this.name in colnames(ubertally)) {
        g <- set.vertex.attribute(g, this.name, value=ubertally[,this.name])
    }

    return(g)

}

## TODO - comment this
## TODO - move the checks to the -test file
OLD.reporting.graph <- function(sim.graph, hidden.popn) {

    sim.settings <- get.graph.attribute(sim.graph, 'sim.settings')
    tau <- sim.settings[['tau']]

    # figure out which vertices are in the hidden popn
    h.idx <- V(sim.graph)[group %in% hidden.popn] 

    ## TODO - LEFT OFF HERE
    ##  - fix bug described below
    ##  - would it be better to make this directed? (a reporting graph is)
    ##  - move checks from comments below to the gnsum-sim-test.Rmd file
    
    # create a reporting graph by remove a fraction of the edges leading to
    # H, according to the true positive rate tau
    pot.edges <- E(sim.graph)[inc(h.idx)]

    ## check: they agree!
    # dbl.edges <- E(sim.graph)[h.idx %--% h.idx]
    # sum(V(sim.graph)[h.idx]$degree)
    # length(pot.edges) + length(dbl.edges)

    # TODO - here is where we need a fix: if we are going to
    #        use tau (and not model reporting some other way),
    #        then strictly speaking we should only remove edges
    #        that go from F to H, not U to H

    # randomly remove edges 
    # NB: as.numeric(pot.edges) convertes the igraph edgelist into a vector
    #     of edge ids, which we then sample
    tolose <- sample(as.numeric(pot.edges), size=(1 - tau)*length(pot.edges))

    rep.graph <- delete.edges(sim.graph, tolose)

    # count up the reports in the reporting graph
    rep.graph <- report.sbm.edges(rep.graph, prefix='y.')    

    ## check: seems to be working!
    #tmp <- get.data.frame(rep.graph, 'vertices')
    ## should be TRUE:
    #all(tmp$y.degree <= tmp$d.degree)
    ## should be FALSE (unless tau==1):
    #all(tmp$y.degree == tmp$d.degree)
    ## y and d should be the same for ties that don't involve H; y should be smaller
    ## for ties that do involve H
    #tmp %>% group_by(group) %>% summarise_each(funs(mean)) %>% select(group, d.FH, y.FH, d.FnotH, y.FnotH)
    
    return(rep.graph)

}

################################################################################
##' 
##' draw.FH.graph
##' 
##' Draw a random graph from a stochastic blockmodel where the groups are
##' determined by membership in the frame population (F) and membership in
##' the hidden population (H). Membership in F is assumed to be
##' independent of membership in H, and vice-versa.
##'
##' @details
##' The four groups, as determined by membership in F and H, are
##' FH, FnotH, notF,notH, notFH.
##'
##' The \code{sim.settings} list should have entries
##' \itemize{
##' \item N the size of the entire population
##' \item p.F the prevalence of the frame population (N.F/N)
##' \item p.H the prevalence of the hidden population (N.H/N)
##' \item pi.within the probability of forming a tie for a pair of
##'       nodes in the same group
##' \item rho.FnotF the relative probability of edge between a node in F
##'       and a node not in F (pi.FnotF = pi.within * rho). so if this is
##'       0.8, then edges between F and not F are 80% as likely as edges between
##'       F and F (or not F and not F)
##' \item rho.HnotH the relative probability of an edge between a node in H
##'       and a node not in H (see rho.FnotF)
##' \item tau the true positive rate
##' }
##'
##' The overall level of connectivity in the graph can be controlled by
##' varying pi.within. The extent to which members of the group are
##' homogenously mixed can be controlled by varying rho.FnotF and rho.HnotH.
##' For example, when rho.FnotF is 1, then the frame population is perfectly
##' mixed with the entire population. When it is less than 1, members of the
##' frame population tend to form ties with one another more than with others.
##' And when it is greater than 1, members of the frame population tend to
##' form ties with others more than with one another.
##'
##' @param sim.settings a list whose entries contain the parameter values
##'        which govern this random graph draw (see details)
##' @return an igraph graph, drawn according to the
##'         settings passed in
draw.FH.graph <- function(sim.settings) {

    # there are more sophisticated ways of doing this, but
    # this has the advantage of being more readable, and
    # causing errors if parameters don't exist
    N <- sim.settings[['N']]
    p.F <- sim.settings[['p.F']]
    p.H <- sim.settings[['p.H']]
    pi.within <- sim.settings[['pi.within']]
    rho.FnotF <- sim.settings[['rho.FnotF']]
    rho.HnotH <- sim.settings[['rho.HnotH']]

    # assume membership in F and H are independent
    N.FnotH <- floor(N * p.F * (1-p.H))
    N.notFnotH <- floor(N * (1-p.F) * (1-p.H))
    N.notFH <- floor(N* (1-p.F) * p.H)
    N.FH <- N - N.FnotH - N.notFnotH - N.notFH

    block.sizes <- c('FnotH' = N.FnotH,
                     'FH' = N.FH,
                     'notFnotH' = N.notFnotH,
                     'notFH' = N.notFH)

    inF <- c(1,1,0,0)
    inH <- c(0,1,0,1)

    # ties between F and not F are either from F to not F or from not F to F:
    Fmask <- (inF %*% t(1-inF)) + ((1-inF) %*% t(inF))
    Fmask <- 1 - Fmask * (1 - rho.FnotF)

    # ties between H and not H are either from H to not H or from not H to H:
    Hmask <- (inH %*% t(1-inH)) + ((1-inH) %*% t(inH))
    Hmask <- 1 - Hmask * (1 - rho.HnotH)

    # first make a matrix of the baseline prob of w/in group edge
    pref.matrix <- matrix(pi.within, nrow=length(block.sizes), ncol=length(block.sizes))
    # then apply the two masks which change probabilities of edges between F / not F and H / not H
    # (assuming these are independent)
    pref.matrix <- pref.matrix * Hmask * Fmask

    this.g <- draw.sbm.graph(block.sizes, pref.matrix)
    this.settings <- sim.settings
    this.res <- set.graph.attribute(this.g, "sim.settings", this.settings)
    return(this.res)

    return(res)
}


