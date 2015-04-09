## ------------------------------------------------------------------------
library(igraph)
library(dplyr)
library(tidyr)

## ------------------------------------------------------------------------
set.seed(12345)

## ------------------------------------------------------------------------
num.sims <- 2

## TODO - parameterize by avg degree w/in group or something?

sim.params <- expand.grid(# size of total population, U
                          N=10e2,
                          #N=10e3,
                          # frame population is half of total
                          p.F=c(1),
                          # hidden popn is 3 percent of total
                          p.H=.5,
                          # prob of edge between two nodes in the same group
                          pi.within=1,
                          # relative prob of edge between two nodes not in the same group
                          rho=c(1),
                          tau=c(.9))



## ------------------------------------------------------------------------
this.g <- draw.4group.graph(sim.params)
this.g.df <- get.data.frame(this.g, 'vertices')

## ------------------------------------------------------------------------
this.r <- reporting.graph(this.g, 
                          hidden.popn=c('FH', 'notFH'),
                          frame.popn=c('FH', 'FnotH'))

## ------------------------------------------------------------------------

these.ests <- reporting.estimates(this.r, 
                                  hidden.popn=c('FH', 'notFH'),
                                  frame.popn=c('FH', 'FnotH'))


## ------------------------------------------------------------------------

system.time(
this.g <- draw.4group.graph(sim.params)
)

system.time(
this.r <- reporting.graph(this.g, 
                          hidden.popn=c('FH', 'notFH'),
                          frame.popn=c('FH', 'FnotH'))
)

this.r.df <- get.data.frame(this.r, 'vertices')

this.r.df %>% group_by(group) %>% summarise(diff=mean(d.degree-y.degree),d.bar=mean(d.degree))

this.r.df %>% group_by(group) %>% summarise(v.bar=mean(v.degree))



