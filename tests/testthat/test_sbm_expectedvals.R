
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

