# =============================================================================
# utils-input.R  - shared input validation for multi-source data
# =============================================================================

# Validate and align a list of data sources measured over the same objects
# (objects in rows). Returns the list with every source reordered to the first
# source's object order, so methods that combine sources positionally (e.g. ADC,
# ADEC) do not silently misalign reordered inputs.
#
# Handles the identifier cases the way a multi-source method needs them:
#   * absent row names  -> alignment is positional; sources must share nrow.
#   * reordered rows     -> reordered to the first source's order.
#   * missing / extra    -> a clear error (these methods need the same objects).
#   * duplicated names   -> a clear error (alignment would be ambiguous).
#
# `what` names the inputs in error messages.
.align_sources <- function(List, what = "the data sources") {
  if (!is.list(List) || length(List) == 0L)
    stop("'List' must be a non-empty list of data matrices.")
  List <- lapply(List, as.matrix)

  nr <- vapply(List, nrow, integer(1))
  if (length(unique(nr)) != 1L)
    stop(sprintf("All of %s must describe the same number of objects (rows); got %s.",
                 what, paste(nr, collapse = ", ")), call. = FALSE)
  if (nr[1] == 0L)
    stop(sprintf("%s are empty (no objects).", what), call. = FALSE)

  rn <- lapply(List, rownames)
  has_names <- vapply(rn, Negate(is.null), logical(1))
  if (!any(has_names)) return(List)            # no identifiers: positional only
  if (!all(has_names))
    stop(sprintf("Either all or none of %s should have row names (object identifiers).",
                 what), call. = FALSE)

  dup <- vapply(rn, anyDuplicated, integer(1))
  if (any(dup > 0L))
    stop(sprintf("%s have duplicated object identifiers (row names).",
                 paste(sprintf("source %d", which(dup > 0L)), collapse = ", ")),
         call. = FALSE)

  ref <- rn[[1]]
  for (k in seq_along(List)) {
    if (setequal(rn[[k]], ref)) {
      if (!identical(rn[[k]], ref))            # same objects, different order
        List[[k]] <- List[[k]][match(ref, rn[[k]]), , drop = FALSE]
    } else {
      miss  <- length(setdiff(ref, rn[[k]]))
      extra <- length(setdiff(rn[[k]], ref))
      stop(sprintf(paste0("source %d does not describe the same objects as source 1 ",
                          "(%d missing, %d extra); %s must share the same objects."),
                   k, miss, extra, what), call. = FALSE)
    }
  }
  List
}
