## For R 2.15.1 and later. Note that calling loadModule() triggers
## a load action, so this does not have to be placed in .onLoad() or evalqOnLoad().

loadModule("dreadds", TRUE)

#' @return A dreaddsModel instance.
dreadds_model <- function(config_file, args = character(0)) {
  return (new(dreaddsModel, config_file, args))
}
