---
title: "Stochastic block model network reporting simulation"
author: "Dennis Feehan"
date: "2015-04-09"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{sbm-simulation}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---




```r
library(igraph)
library(dplyr)
library(tidyr)
library(stringr)
```


```r
set.seed(12345)
```

Now set up the params for one random draw.


```r
num.sims <- 2

## TODO - parameterize by avg degree w/in group or something?

sim.params <- expand.grid(# size of total population, U
                          N=10e2,
                          #N=10e3,
                          # frame population is half of total
                          p.F=1,
                          # hidden popn is 3 percent of total
                          p.H=.5,
                          # prob of edge between two nodes in the same group
                          pi.within=.05,
                          # relative prob of edge between two nodes not in the same group
                          rho=c(1),
                          tau=c(.9))
```

Now draw the actual random graph.

TODO - LEFT OFF HERE

* converting sim code to package
* I think this test file has slightly old version of code
* code up sim using the package (specify version)
* re-run sim to specify according to writeup, if necessary

TODO - TO DEMONSTRATE

* report.sbm.edges

TODO - UNIT TESTS

* report.sbm.edges (be sure to mention counter-intuitive in vs out arg)
* reporting.estimates

etc...


```r
this.g <- draw.4group.graph(sim.params, type="nested")
this.g.df <- get.data.frame(this.g, 'vertices')
```

Next, generate the reporting graph.


```r
this.r <- reporting.graph(this.g, 
                          hidden.popn=c('FH', 'notFH'),
                          frame.popn=c('FH', 'FnotH'))
```

Finally, produce estimates


```r
these.ests <- reporting.estimates(this.r, 
                                  hidden.popn=c('FH', 'notFH'),
                                  frame.popn=c('FH', 'FnotH'))
```


```r
#system.time(
#this.g <- draw.4group.graph(sim.params)
#)

#system.time(
#this.r <- reporting.graph(this.g, 
#                          hidden.popn=c('FH', 'notFH'),
#                          frame.popn=c('FH', 'FnotH'))
#)

this.r.df <- get.data.frame(this.r, 'vertices')

this.r.df %>% group_by(group) %>% summarise(diff=mean(d.degree-y.degree),d.bar=mean(d.degree))
```

```
## Source: local data frame [2 x 3]
## 
##   group  diff  d.bar
## 1    FH 5.096 50.974
## 2 FnotH 0.000 49.946
```

```r
this.r.df %>% group_by(group) %>% summarise(v.bar=mean(v.degree))
```

```
## Source: local data frame [2 x 2]
## 
##   group  v.bar
## 1    FH 48.292
## 2 FnotH 47.532
```

## Code checks

TODO - these should be made into unit tests...

Here are a few internal consistency checks that I ran to be sure the
simulation code is working.

This should be TRUE:

```r
this.g.df <- get.data.frame(this.g, 'vertices')

# if p.F = 1, then this 
if (! "d.notFnotH" %in% colnames(this.g.df)) {
    this.g.df$d.notFnotH <- 0
}

with(this.g.df, all(d.degree == d.FnotH + d.FH + d.notFnotH + d.notFH))
```

```
## [1] TRUE
```

```r
with(this.g.df, all(d.degree == d.FnotH + d.FH + d.notFH))
```

```
## [1] TRUE
```

The `diff` column of the dataframe produced here should all be zeroes.










