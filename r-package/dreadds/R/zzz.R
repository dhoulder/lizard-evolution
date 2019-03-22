## For R 2.15.1 and later. Note that calling loadModule() triggers
## a load action, so this does not have to be placed in .onLoad() or evalqOnLoad().

loadModule("dreadds", TRUE)

.build_args <- function(args, ...) {
    for (ov in list(...)) {
        if (length(ov) < 2)
            next
        args <- c(args, ov)
    }
    return (args)
}


#' @return A dreaddsModel instance.
createDreadDS <- function(config.file = NULL,
                          genetic.dims = NULL,
                          gene.flow.max.distance = NULL,
                          niche.evolution.rate = NULL,
                          dispersal.min = NULL,
                          output.dir = NULL,
                          output.file.prefix = NULL,
                          iterations = NULL,
                          check.speciation = NULL,
                          args = character(0)) {

    av <- .build_args(args,
                      c('-c', config.file),
                      c('--genetic-dims',  genetic.dims),
                      c('--gene-flow-max-distance',  gene.flow.max.distance),
                      c('--niche-evolution-rate', niche.evolution.rate),
                      c('--dispersal-min', dispersal.min),
                      c('-o', output.dir),
                      c('--output-file-prefix', output.file.prefix),
                      c('-n', iterations),
                      c('--check-speciation', check.speciation))
    return (new(dreaddsModel, av))
}
