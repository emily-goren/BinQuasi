#'
#' Extract model coefficients
#' 
#' @description Extract model coefficients.
#' 
#' @param object The fitted model.
#' @param type Must be one of \code{raw}, \code{bias}, or \code{corrected}.
#' @param ... Additional arguments.
#' 
#' @rdname coef.glm
#' @export
#'
coef.glm=function(object, type=c('raw','bias','corrected'), ...)
{
  type=match.arg(type)
  if(type=='raw') return(object$coefficients)
  
  good.wt=object$weights>0
  this.qr=object$qr
  this.hatd=.rowSums(qr.Q(this.qr)[,seq_len(this.qr$rank), drop=FALSE]^2, sum(good.wt), this.qr$rank)## not affected by pivoting
  
  this.mu=object$fitted.values[good.wt]
  this.eta=object$family$linkfun(this.mu)
  if (is.null(object$family$d2linkfun))
    object$family$d2linkfun <- function(mu) -1/mu/mu
  this.w2ksi= 0.5* this.hatd / sqrt(object$weights[good.wt]) *object$family$d2linkfun(this.mu)*object$family$mu.eta(this.eta)^2  ## this requires existence of d2linkfun components on the family
  this.bias=qr.coef(this.qr, this.w2ksi)
  if(type=='bias') return(this.bias)
  
  ans=object$coefficient-drop(this.bias)
  attr(ans, 'bias')=this.bias
  ans
}
