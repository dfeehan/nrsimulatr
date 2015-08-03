
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


