## -----------------------------------------------------------------------------
library(knitr)
opts_chunk$set(eval=FALSE)

## -----------------------------------------------------------------------------
#  library(igraph)
#  library(dplyr)
#  library(tidyr)
#  library(stringr)
#  
#  library(nrsimulatr)

## -----------------------------------------------------------------------------
#  set.seed(12345)

## -----------------------------------------------------------------------------
#  sim.params <- expand.grid(# size of total population, U
#                            N=10e3,
#                            # fraction of the population that is in the frame population (F)
#                            p.F=.5,
#                            # fraction of the population that is in the hidden popn (H)
#                            p.H=.03,
#                            # probability of being in F, conditional on being in H
#                            p.F.given.H=1,
#                            # prob of edge between two nodes in the same group
#                            zeta=.05,
#                            # relative prob of edge between hidden and not hidden
#                            rho=c(.5, .8, 1),
#                            # relative prob of edge between frame and not frame
#                            xi = 0.6,
#                            # sampling fraction for simple random sample w/out replacement from F
#                            frame.sampling.frac=.3,
#                            # sampling fraction for relative probability sample from H
#                            hidden.sampling.frac=.4,
#                            # true positive rate
#                            tau=.8)
#  
#  # add an id for each parameter combination
#  sim.params$param.index <- 1:nrow(sim.params)
#  

## -----------------------------------------------------------------------------
#  
#  all.param.list <- split(sim.params, rownames(sim.params))
#  
#  # we would change 'type' if we wanted to use a different type of
#  # simulation model -- for example, a one-parameter stochastic blockmodel with four groups
#  all.param.list <- plyr::llply(all.param.list,
#                                sbm_params,
#                                type="4group_2param")
#  

## -----------------------------------------------------------------------------
#  
#  ## we'll repeat each set of parameter combinations this number of times
#  ## (small here to keep this reasonably fast; in practice, you probably want many more)
#  num.sims <- 2
#  
#  sim.ests <- plyr::ldply(all.param.list,
#               function(these.sim.params) {
#  
#                  # construct an object containing the parameters that govern
#                  # reporting (in our case, this will use our tau parameter to
#                  # simulate false negatives)
#                  rep.params <- imperfect_reporting(these.sim.params)
#  
#                  # construct an object containing the parameters that govern
#                  # sampling from the frame population
#                  frame.sample.params <- frame_srs(list(sampling.frac=
#                                                        these.sim.params$frame.sampling.frac))
#  
#                  # construct an object containing the parameters that govern
#                  # sampling from the hidden population
#                  hidden.sample.params <- hidden_relprob(list(sampling.frac=
#                                                              these.sim.params$hidden.sampling.frac,
#                                                              inclusion.trait='d.degree'))
#  
#                  this.N <- these.sim.params$N
#  
#                  ## this goes through each repetition for the parameter set we're on
#                  these.res <- plyr::ldply(1:num.sims,
#                                     function(this.idx) {
#  
#                                       # draw a random graph from the stochastic block model
#                                       this.g <- generate_graph(these.sim.params)
#  
#                                       # use the reporting parameters to convert the
#                                       # undirected social network into a directed
#                                       # reporting graph
#                                       this.r <- reporting_graph(rep.params, this.g)
#  
#                                       # get some quantities about the entire reporting graph
#                                       # from a census
#                                       census.df <- sample_graph(entire_census(), this.r)
#                                       #res <- sbm_census_quantities(census.df)
#  
#                                       res <- list()
#  
#                                       # get sample-based estimates from
#                                       # a sample of the frame population and a sample
#                                       # of the hidden population
#                                       frame.df <- sample_graph(frame.sample.params,
#                                                                this.r)
#                                       hidden.df <- sample_graph(hidden.sample.params,
#                                                                 this.r)
#  
#                                       # note that we will assume that the degrees
#                                       # are correct (ie, we are not simulating the
#                                       # known population method)
#                                       res$sample.y.F.H <- with(frame.df,
#                                                                sum((y.FH + y.notFH)*
#                                                                    sampling.weight))
#  
#                                       res$sample.d.F.U <- with(frame.df,
#                                                                sum(d.degree*sampling.weight))
#  
#                                       res$sample.dbar.F.F <- with(frame.df,
#                                                                   sum((d.FH + d.FnotH)*
#                                                                       sampling.weight)/
#                                                                   sum(sampling.weight))
#  
#                                       res$sample.vbar.H.F <- with(hidden.df,
#                                                                sum((v.FnotH + v.FH)*
#                                                                    sampling.weight)/
#                                                                sum(sampling.weight))
#  
#                                       res$sample.basic <- with(res,
#                                                                (sample.y.F.H / sample.d.F.U)*
#                                                                this.N)
#  
#                                       res$sample.adapted <- with(res,
#                                                                  sample.y.F.H / sample.dbar.F.F)
#  
#                                       res$sample.generalized <- with(res,
#                                                                      sample.y.F.H / sample.vbar.H.F)
#  
#                                       # for convenience, keep track of which rep this is
#                                       res$rep <- this.idx
#  
#                                       return(data.frame(res))
#                                     })
#  
#                  ## for convenience, keep track of which parameter set this is
#                  these.res$param.index <- these.sim.params$param.index
#  
#                  return(these.res)
#               })
#  

## -----------------------------------------------------------------------------
#  
#  ci_low <- function(x) { quantile(x, .025, na.rm=TRUE) }
#  ci_high <- function(x) { quantile(x, .975, na.rm=TRUE) }
#  
#  sim.agg.byest <- sim.ests %>%
#                   group_by(param.index) %>%
#                   rename(basic=sample.basic,
#                          adapted=sample.adapted,
#                          generalized=sample.generalized) %>%
#                   summarise_each(funs(mean, ci_low, ci_high),
#                                  basic, adapted, generalized) %>%
#                   gather(qty, value,
#                          matches("basic|adapted|generalized")) %>%
#                   separate(qty, c("estimator", "qty"), sep="\\_", extra="merge") %>%
#                   spread(qty, value)
#  
#  # merge in parameter values
#  sim.params <- sim.params %>% mutate(N.H = p.H * N)
#  
#  sim.agg.byest <- inner_join(sim.agg.byest, sim.params, by='param.index')
#  

## ---- results='asis'----------------------------------------------------------
#  
#  sim.agg.byest <- sim.agg.byest %>% mutate(bias = mean - N.H)
#  
#  bias.tab <- sim.agg.byest %>% select(estimator, rho, bias) %>%
#                                spread(estimator, bias)
#  
#  kable(bias.tab)
#  

