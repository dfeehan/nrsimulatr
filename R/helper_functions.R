
#' assign labels to vertices in a graph
#'
#' @param graph the \code{igraph} graph
#' @param block.sizes a vector with the number of nodes in each group
#' @param block.names a vector (same length as \code{block.sizes}) with names of each group
#' @export
label.vertices <- function(graph, block.sizes, block.names) {

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

