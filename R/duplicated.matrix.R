#'
#' Determine duplicate elements
#' 
#' @description Determine which elements are duplicates.
#' 
#' @param x A matrix.
#' @param incomparables The values that cannot be compared. FALSE means all values can be compared.
#' @param MARGIN The array margin to hold fixed.
#' @param fromLast Logical. Should duplication be considered from teh last instance?
#' @param ... Additional arguments.
#' 
#' @export
#'
duplicated.matrix = function(x, incomparables = FALSE, MARGIN = 1L, fromLast = FALSE, ...)
{
  if (!is.matrix(x) ||!is.numeric(x[1]) || !identical(incomparables, FALSE) || MARGIN!=1L)
    return(base::duplicated.matrix(x, incomparables, MARGIN, fromLast, ...))
  .Call(C_dupRowNumMat,x, as.logical(fromLast))
}

#'
#' Remove duplicate elements
#' 
#' @description Remove duplicate elements.
#' 
#' @param x A matrix.
#' @param incomparables The values that cannot be compared. FALSE means all values can be compared.
#' @param MARGIN The array margin to hold fixed.
#' @param fromLast Logical. Should duplication be considered from teh last instance?
#' @param ... Additional arguments.
#' 
#' @export
#' 
unique.matrix = function(x, incomparables = FALSE, MARGIN = 1, fromLast = FALSE, ...)
{
  if (!is.matrix(x) ||!is.numeric(x[1]) || !identical(incomparables, FALSE) || MARGIN!=1L)
    return(base::unique.matrix(x, incomparables, MARGIN, fromLast, ...))
  x[!.Call(C_dupRowNumMat,x,as.logical(fromLast)),,drop=FALSE]
}
