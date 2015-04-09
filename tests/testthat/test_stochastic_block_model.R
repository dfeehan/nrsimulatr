
context("sbm preference matrix")

test.block.sizes <- c("FnotH"=100, "FH"=100, "notFnotH"=100, "notFH"=100)
test.inF <- c(1, 1, 0, 0)
test.inH <- c(0, 1, 0, 1)

##########
## matrix with all entries the same (simple)
res <- pref.matrix.4group(block.sizes=test.block.sizes,
                          pi.within=0.05, rho=1, type="simple")

expect_equal(dim(res), c(4,4))
expect_true(all(res==0.05))

##########
## matrix with all entries the same (nested)
res <- pref.matrix.4group(block.sizes=test.block.sizes,
                          pi.within=0.05, rho=1, type="nested",
                          inF=test.inF, inH=test.inH)

expect_equal(dim(res), c(4,4))
expect_true(all(res==0.05))

##########
## inhomogenous mixing (simple)

res <- pref.matrix.4group(block.sizes=test.block.sizes,
                          pi.within=0.05, rho=.1, type="simple")

expect_equal(dim(res), c(4,4))
expect_true(all(diag(res)==0.05))

expect_equal(setNames(res[,1],NULL),
             c(0.05, 0.005, 0.005, 0.005),
             tolerance=.001)

##########
## inhomogenous mixing (simple, blocks reordered)

res <- pref.matrix.4group(block.sizes=test.block.sizes[c(3,1,4,2)],
                          pi.within=0.05, rho=.1, type="simple")

expect_equal(dim(res), c(4,4))
expect_true(all(diag(res)==0.05))

expect_equal(setNames(res[,1],NULL),
             c(0.05, 0.005, 0.005, 0.005),
             tolerance=.001)

##########
## inhomogenous mixing (nested)

res <- pref.matrix.4group(block.sizes=test.block.sizes,
                          pi.within=0.05, rho=.1, type="nested",
                          inF=test.inF, inH=test.inH)

expect_equal(dim(res), c(4,4))
expect_true(all(diag(res)==0.05))

expect_equal(setNames(res[,1],NULL),
             c(0.05, 0.005, 0.005, 0.0005),
             tolerance=.001)

## check that the reverse diagonal (bottom-left to top-right)
## is all rho^2
expect_equal(diag(apply(res, 1, rev)),
             c(.0005, .0005, .0005, .0005))


##########
## inhomogenous mixing (nested, blocks reordered)

res <- pref.matrix.4group(block.sizes=test.block.sizes[c(3,1,4,2)],
                          pi.within=0.05, rho=.1, type="nested",
                          inF=test.inF[c(3,1,4,2)], 
                          inH=test.inH[c(3,1,4,2)])

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

test.block.sizes <- c("FnotH"=100, "FH"=100, "notFnotH"=100, "notFH"=100)
test.inF <- c(1, 1, 0, 0)
test.inH <- c(0, 1, 0, 1)
test.inFnotH <- c(1,0,0,0)

##########
## matrix with all entries the same (simple)
test.pref.matrix <- pref.matrix.4group(block.sizes=test.block.sizes,
                                       pi.within=0.05, rho=1, type="simple")

## edges within a block
eec <- expected.edgecount.4group(test.pref.matrix,
                                 block.sizes=test.block.sizes,
                                 first.group=c(1,0,0,0),
                                 second.group=c(1,0,0,0))

expect_equal(eec, choose(100,2)*.05)


## edges between two distinct blocks
eec <- expected.edgecount.4group(test.pref.matrix,
                                 block.sizes=test.block.sizes,
                                 first.group=c(1,0,0,0),
                                 second.group=c(0,1,0,0))

expect_equal(eec, 100*100*.05)

## edges between two overlapping blocks
eec <- expected.edgecount.4group(test.pref.matrix,
                                 block.sizes=test.block.sizes,
                                 first.group=c(1,1,0,0),
                                 second.group=c(0,1,0,0))

expect_equal(eec, 100*100*.05 + choose(100,2)*.05)

## edges between two overlapping blocks
eec <- expected.edgecount.4group(test.pref.matrix,
                                 block.sizes=test.block.sizes,
                                 first.group=c(0,1,0,0),
                                 second.group=c(1,1,0,0))

expect_equal(eec, 100*100*.05 + choose(100,2)*.05)

## edges between two overlapping blocks
eec <- expected.edgecount.4group(test.pref.matrix,
                                 block.sizes=test.block.sizes,
                                 first.group=c(1,1,0,0),
                                 second.group=c(1,1,0,0))

expect_equal(eec, 2*(100*100*.05 + choose(100,2)*.05))

##########
## matrix with inhomogenous mixing (simple)
test.pref.matrix <- pref.matrix.4group(block.sizes=test.block.sizes,
                                       pi.within=0.05, rho=.1, type="simple")

## edges between two overlapping blocks
eec <- expected.edgecount.4group(test.pref.matrix,
                                 block.sizes=test.block.sizes,
                                 first.group=c(1,1,0,0),
                                 second.group=c(1,1,0,0))

expect_equal(eec, 2*(100*100*.05*.1 + choose(100,2)*.05))

##########
## matrix with inhomogenous mixing (nested)
test.pref.matrix <- pref.matrix.4group(block.sizes=test.block.sizes,
                                       pi.within=0.05, rho=.1, 
                                       type="nested",
                                       inF=test.inF, inH=test.inH)

## edges between two overlapping blocks
eec <- expected.edgecount.4group(test.pref.matrix,
                                 block.sizes=test.block.sizes,
                                 first.group=c(1,0,0,1),
                                 second.group=c(0,0,0,1))

expect_equal(eec, 100*100*.05*.1*.1 + choose(100,2)*.05)


########################################
context("sbm draw graph")

test.sim.params <- list(N=10000,
                        p.F=.5,
                        p.H=.1,
                        pi.within=.1,
                        rho=.5,
                        tau=1)

test.g <- draw.4group.graph(sim.settings = test.sim.params)

N.FnotH <- with(test.sim.params, floor(N * p.F * (1-p.H)))
N.notFnotH <- with(test.sim.params, floor(N * (1-p.F) * (1-p.H)))
N.notFH <- with(test.sim.params, floor(N * (1-p.F) * p.H))
N.FH <- with(test.sim.params, N - N.FnotH - N.notFnotH - N.notFH)

test.block.sizes <- c('FnotH' = N.FnotH,
                      'FH' = N.FH,
                      'notFnotH' = N.notFnotH,
                      'notFH' = N.notFH)

test.pref.matrix <- with(test.sim.params,
                         pref.matrix.4group(block.sizes=test.block.sizes,
                                            pi.within=pi.within,
                                            rho=rho,
                                            type="simple"))

test.F.df <- get.data.frame(test.g, "vertices") %>% filter(group %in% c("FnotH", "FH")) 

test.res <- test.F.df %>% 
            summarise(deg.tot = sum(d.degree),
                      d.bar.F.FnotH = mean(d.FnotH),
                      d.F.FH = sum(d.FH),
                      d.bar.F.F = mean(d.FnotH + d.FH),
                      d.F.notFnotH = sum(d.notFnotH))

test.inF <- c(1,1,0,0)

## NB: there will be two reports for each edge...
eec.d.bar.F.FnotH <- expected.edgecount.4group(test.pref.matrix,
                                               block.sizes=test.block.sizes,
                                               first.group=test.inF,
                                               second.group=c(1,0,0,0)) /
                     sum(test.block.sizes*test.inF)

eec.d.F.F <- expected.edgecount.4group(test.pref.matrix,
                                       block.sizes=test.block.sizes,
                                       first.group=test.inF,
                                       second.group=test.inF)

## TODO - LEFT OFF HERE -- finish this set of test (for random network draws)



