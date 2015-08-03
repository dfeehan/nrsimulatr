
########################################
context("sbm reporting edges")

#g <- graph.formula( Ace-Bob-Cecil-Alice, Daniel-Cecil-Eugene, Cecil-Gordon )
g <- graph.formula( A+-+B-+C+-+A, D+-C+-+E, C+-G )

g <- set.vertex.attribute(g,
                          "group",
                          value=get.vertex.attribute(g, "name"))

g <- set.graph.attribute(g,
                         "block.sizes",
                         value=c("A"=1, "B"=1, "C"=1, "D"=1, "E"=1, "F"=1, "G"=1))

# hack this into an sbm object to test the report_edges.sbm function
class(g) <- append(class(g), "sbm")

# you can run this to visualize this graph if you want to check these tests
# plot(g)

# for now, mode='in' still means out-reports, but this will likely change
res.outrep <- report_edges(g, prefix='y.', mode='in')
rep.df <- get.data.frame(res.outrep, 'vertices')

test_that("counting out-reports is accurate", {
    expect_equal(rep.df$y.degree, c(2,1,4,1,1,0))
    expect_equal(rep.df$y.A, c(0,1,1,0,0,0))
    expect_equal(rep.df$y.B, c(1,0,0,0,0,0))
    expect_equal(rep.df$y.C, c(1,1,0,0,1,1))
    expect_equal(rep.df$y.D, c(0,0,1,0,0,0))
    expect_equal(rep.df$y.E, c(0,0,1,0,0,0))
    expect_equal(rep.df$y.F, c(0,0,0,0,0,0))
    expect_equal(rep.df$y.G, c(0,0,0,0,0,0))
})

res.vis <- report_edges(g, prefix='v.', mode='out')
vis.df <- get.data.frame(res.vis, 'vertices')

test_that("counting in-reports is accurate", {
    expect_equal(vis.df$v.degree, c(2,2,3,0,1,1))
    expect_equal(vis.df$v.A, c(0,1,1,0,0,0))
    expect_equal(vis.df$v.B, c(1,0,1,0,0,0))
    expect_equal(vis.df$v.C, c(1,0,0,1,1,0))
    expect_equal(vis.df$v.D, c(0,0,0,0,0,0))
    expect_equal(vis.df$v.E, c(0,0,1,0,0,0))
    expect_equal(vis.df$v.F, c(0,0,0,0,0,0))
    expect_equal(vis.df$v.G, c(0,0,1,0,0,0))
})

## TODO - more tests could be sure aggregates work
## (i.e., group vertices above into two groups and check again)



