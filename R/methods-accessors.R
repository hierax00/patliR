#' @include AllGenerics.R
NULL

#' @rdname compounds
setMethod("compounds", "PatliRProject", function(proj) proj@compounds)

#' @rdname compounds
setMethod("compounds<-", "PatliRProject", function(proj, value) {
  proj@compounds <- value
  methods::validObject(proj)
  proj
})

#' @rdname matrixRaw
setMethod("matrixRaw", "PatliRProject", function(proj) proj@matrix_raw)

#' @rdname matrixRaw
setMethod("matrixRaw<-", "PatliRProject", function(proj, value) {
  proj@matrix_raw <- value
  methods::validObject(proj)
  proj
})

#' @rdname binarizedMatrix
setMethod("binarizedMatrix", "PatliRProject", function(proj) proj@binarized)

#' @rdname binarizedMatrix
setMethod("binarizedMatrix<-", "PatliRProject", function(proj, value) {
  proj@binarized <- value
  methods::validObject(proj)
  proj
})

#' @rdname projectLog
setMethod("projectLog", "PatliRProject", function(proj) proj@log)

#' @rdname projectLog
setMethod("projectLog<-", "PatliRProject", function(proj, value) {
  proj@log <- value
  methods::validObject(proj)
  proj
})

#' @rdname projectDir
setMethod("projectDir", "PatliRProject", function(proj) proj@project_dir)

#' @rdname projectDir
setMethod("cacheDir", "PatliRProject", function(proj) proj@cache_dir)

#' @rdname patliRResults
setMethod("patliRResults", "PatliRProject", function(proj, name = NULL) {
  if (is.null(name)) return(proj@results)
  proj@results[[name]]
})

#' @rdname patliRResults
setMethod("patliRResults<-", "PatliRProject", function(proj, name, value) {
  stopifnot(is.character(name), length(name) == 1, nzchar(name))
  proj@results[[name]] <- value
  methods::validObject(proj)
  proj
})
