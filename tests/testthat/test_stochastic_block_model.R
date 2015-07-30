
context("sbm preference matrix")

##########
## matrix with all entries the same (simple)
sim.params <- list(N=400, p.F=.5, p.H=.5, pi.within=.05, rho=1)
params.simple <- sbm_params(sim.params, type="4group_1param_simple")

res <- sbm_pref_matrix(params.simple)

expect_equal(dim(res), c(4,4))
expect_true(all(res==0.05))

##########
## matrix with all entries the same (nested)
sim.params <- list(N=400, p.F=.5, p.H=.5, pi.within=.05, rho=1)
params.nested <- sbm_params(sim.params, type="4group_1param_nested")

res <- sbm_pref_matrix(params.nested)

expect_equal(dim(res), c(4,4))
expect_true(all(res==0.05))

##########
## inhomogenous mixing (simple)
sim.params <- list(N=400, p.F=.5, p.H=.5, pi.within=.05, rho=.1)

params.simple <- sbm_params(sim.params, type="4group_1param_simple")

res <- params.simple$pref.matrix

expect_equal(dim(res), c(4,4))
expect_true(all(diag(res)==0.05))

expect_equal(setNames(res[,1],NULL),
             c(0.05, 0.005, 0.005, 0.005),
             tolerance=.001)

##########
## inhomogenous mixing (nested)
sim.params <- list(N=400, p.F=.5, p.H=.5, pi.within=.05, rho=.1)

params.nested <- sbm_params(sim.params, type="4group_1param_nested")

res <- params.nested$pref.matrix

expect_equal(dim(res), c(4,4))
expect_true(all(diag(res)==0.05))

expect_equal(setNames(res[,1],NULL),
             c(0.05, 0.005, 0.005, 0.0005),
             tolerance=.001)

## check that the reverse diagonal (bottom-left to top-right)
## is all rho^2
expect_equal(diag(apply(res, 1, rev)),
             c(.0005, .0005, .0005, .0005))


########################################
context("sbm expected number of edges")

##########
## matrix with all entries the same (simple)
sim.params <- list(N=400, p.F=.5, p.H=.5, pi.within=.05, rho=1)
params.simple <- sbm_params(sim.params, type="4group_1param_simple")

# order: 'FnotH', 'FH', 'notFnotH', 'notFH'

## edges within a block
eec <- expected_edgecount(params.simple,
                          first.group='FnotH',
                          second.group='FnotH')

expect_equal(eec, choose(100,2)*.05)


## edges between two distinct blocks
eec <- expected_edgecount(params.simple,
                          first.group='FnotH',
                          second.group='FH')

expect_equal(eec, 100*100*.05)

## edges between two overlapping blocks
eec <- expected_edgecount(params.simple,
                          first.group=c('FnotH', 'FH'),
                          second.group='FH')

expect_equal(eec, 100*100*.05 + choose(100,2)*.05)

## edges between two overlapping blocks
eec <- expected_edgecount(params.simple,
                          first.group='FH',
                          second.group=c('FnotH', 'FH'))

expect_equal(eec, 100*100*.05 + choose(100,2)*.05)

## edges between two overlapping blocks
eec <- expected_edgecount(params.simple,
                          first.group=c('FnotH', 'FH'),
                          second.group=c('FnotH', 'FH'))
expect_equal(eec, 2*(100*100*.05 + choose(100,2)*.05))

##########
## matrix with inhomogenous mixing (simple)
sim.params <- list(N=400, p.F=.5, p.H=.5, pi.within=.05, rho=.1)
params.simple <- sbm_params(sim.params, type="4group_1param_simple")

## edges between two overlapping blocks
eec <- expected_edgecount(params.simple,
                          first.group=c('FnotH', 'FH'),
                          second.group=c('FnotH', 'FH'))

expect_equal(eec, 2*(100*100*.05*.1 + choose(100,2)*.05))

##########
## matrix with inhomogenous mixing (nested)
sim.params <- list(N=400, p.F=.5, p.H=.5, pi.within=.05, rho=.1)
params.nested <- sbm_params(sim.params, type="4group_1param_nested")

## edges between two overlapping blocks
eec <- expected_edgecount(params.nested,
                          first.group=c('FnotH', 'notFH'),
                          second.group='notFH')

expect_equal(eec, 100*100*.05*.1*.1 + choose(100,2)*.05)

########################################
context("sbm draw graph: 4group_1param_simple")

# given an sbm graph, 
# count the number of connections between two given blocks
get.conns <- function(g, gp1, gp2) {

    gp1.idx <- V(g)$group == gp1
    gp2.idx <- V(g)$group == gp2
    
    tot <- length(E(g)[gp1.idx %--% gp2.idx])
    
    return(data.frame(from=gp1, to=gp2,
                      tot=tot,
                      mean=tot/length(gp1.idx)))

}

# get the expected number of connections between blocks based
# on parameter values
get.block.evs <- function(params) {
    # get a table of expected values
    gps <- names(params$block.sizes)
    combs <- expand.grid(from=gps, to=gps)

    res <- adply(combs,
                 1,
                 function(x) {
                     this.ev <- expected_edgecount(params, x$from, x$to)
                     return(data.frame(from=x$from, to=x$to, expected.tot=this.ev))
                 })

    return(res)
}

# count the actual connections between blocks 
# in a randomly sampled graph
get.block.conns <- function(g, params) {
    # get a table of expected values
    gps <- names(params$block.sizes)
    combs <- expand.grid(from=gps, to=gps)

    res <- adply(combs,
                 1,
                 function(x) {
                     this.conns <- get.conns(g, x$from, x$to)
                     return(this.conns)
                 })

    return(res)
}


##### now start the actual testing

## 4group_1param_simple

test.sim.params <- list(N=500,
                        p.F=.5,
                        p.H=.3,
                        pi.within=.1,
                        xi=.5,
                        rho=.5,
                        tau=1)

params <- sbm_params(test.sim.params, type="4group_1param_simple")

test.ev <- get.block.evs(params)

M <- 100
res <- ldply(1:M,
             function(x) {
                 test.g <- generate_graph(params)
                 return(get.block.conns(test.g,params))
             },
             .progress='text')

# compare the M sampled graphs to the expected values
res.summ <- res %>% group_by(from,to) %>% summarise(mean.tot = mean(tot),
                                                    sd.tot = sd(tot))

res.both <- join(res.summ, test.ev)

test_that("4group_1param_simple graphs get generated correctly", {
   # we use tolerance=.1 because we're comparing simulated data to
   # expected values, and we wouldn't expect precise agreement
   expect_equal(res.both$mean.tot, res.both$expected.tot, tolerance=.1)
})

## 4group_1param_nested

test.sim.params <- list(N=500,
                        p.F=.5,
                        p.H=.3,
                        pi.within=.1,
                        xi=.5,
                        rho=.5,
                        tau=1)

params <- sbm_params(test.sim.params, type="4group_1param_nested")

test.ev <- get.block.evs(params)

M <- 100
res <- ldply(1:M,
             function(x) {
                 test.g <- generate_graph(params)
                 return(get.block.conns(test.g,params))
             },
             .progress='text')

# compare the M sampled graphs to the expected values
res.summ <- res %>% group_by(from,to) %>% summarise(mean.tot = mean(tot),
                                                    sd.tot = sd(tot))

res.both <- join(res.summ, test.ev)

test_that("4group_1param_nested graphs get generated correctly", {
   # we use tolerance=.1 because we're comparing simulated data to
   # expected values, and we wouldn't expect precise agreement
   expect_equal(res.both$mean.tot, res.both$expected.tot, tolerance=.1)
})

## 4group_2param

test.sim.params <- list(N=500,
                        p.F=.5,
                        p.H=.3,
                        pi.within=.1,
                        xi=.5,
                        rho=.5,
                        tau=1)

params <- sbm_params(test.sim.params, type="4group_2param")

test.ev <- get.block.evs(params)

M <- 100
res <- ldply(1:M,
             function(x) {
                 test.g <- generate_graph(params)
                 return(get.block.conns(test.g,params))
             },
             .progress='text')

# compare the M sampled graphs to the expected values
res.summ <- res %>% group_by(from,to) %>% summarise(mean.tot = mean(tot),
                                                    sd.tot = sd(tot))

res.both <- join(res.summ, test.ev)

test_that("4group_2param", {
   # we use tolerance=.1 because we're comparing simulated data to
   # expected values, and we wouldn't expect precise agreement
   expect_equal(res.both$mean.tot, res.both$expected.tot, tolerance=.1)
})

########################################
context("sbm reporting edges")

# TODO

