#' @include AllClasses.R
NULL

#' Access the compounds table of a patliR project
#'
#' @param proj A [`PatliRProject-class`] object.
#' @return A `data.frame` with columns `id`, `pubchem_id`, `name`, `smiles`,
#'   `canonical_smiles`, `source`; zero rows if [prep_compounds()] has not
#'   run yet.
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compounds(proj)
#' @export
setGeneric("compounds", function(proj) standardGeneric("compounds"))

#' @rdname compounds
#' @param value Replacement `data.frame`, validated on assignment.
#' @export
setGeneric("compounds<-", function(proj, value) standardGeneric("compounds<-"))

#' Access the raw (non-binarized) abundance matrix of a patliR project
#'
#' @inheritParams compounds
#' @return A `data.frame` with one row per compound and one column per
#'   condition, or an empty `data.frame` if [prep_binarize()] has not run.
#' @export
setGeneric("matrixRaw", function(proj) standardGeneric("matrixRaw"))

#' @rdname matrixRaw
#' @inheritParams compounds<-
#' @export
setGeneric("matrixRaw<-", function(proj, value) standardGeneric("matrixRaw<-"))

#' Access the binarized (presence/absence) matrix of a patliR project
#'
#' @inheritParams compounds
#' @return A `data.frame` with the same shape as [matrixRaw()], with 0/1
#'   calls, or an empty `data.frame` if [prep_binarize()] has not run.
#' @export
setGeneric("binarizedMatrix", function(proj) standardGeneric("binarizedMatrix"))

#' @rdname binarizedMatrix
#' @inheritParams compounds<-
#' @export
setGeneric("binarizedMatrix<-", function(proj, value) standardGeneric("binarizedMatrix<-"))

#' Access the accumulated event log of a patliR project
#'
#' @inheritParams compounds
#' @return A `data.frame` with columns `step`, `id`, `message`, `timestamp`.
#' @export
setGeneric("projectLog", function(proj) standardGeneric("projectLog"))

#' @rdname projectLog
#' @inheritParams compounds<-
#' @export
setGeneric("projectLog<-", function(proj, value) standardGeneric("projectLog<-"))

#' Access the project and cache directories of a patliR project
#'
#' @inheritParams compounds
#' @return A character scalar with the absolute path.
#' @export
setGeneric("projectDir", function(proj) standardGeneric("projectDir"))

#' @rdname projectDir
#' @export
setGeneric("cacheDir", function(proj) standardGeneric("cacheDir"))

#' Access one entry of the flexible results bag of a patliR project
#'
#' @description
#' `patliRResults()` is the read accessor for everything that is not one of
#' the structural slots of [`PatliRProject-class`] (i.e. everything besides
#' `compounds`, `matrix_raw`, `binarized`): ADME tables, target predictions,
#' one network per condition, bias audits, rankings, and so on. Each family
#' of functions documents which name(s) it writes to this bag.
#'
#' @inheritParams compounds
#' @param name `NULL` (the default) to get the whole named list, or a
#'   character scalar to get one entry (`NULL` if that entry does not exist
#'   yet).
#' @return A named `list`, or one element of it.
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' patliRResults(proj) # empty list: nothing has been computed yet
#' @export
setGeneric("patliRResults", function(proj, name = NULL) standardGeneric("patliRResults"))

#' @rdname patliRResults
#' @param value Replacement value for the entry named `name`.
#' @export
setGeneric("patliRResults<-", function(proj, name, value) standardGeneric("patliRResults<-"))
