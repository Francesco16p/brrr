#' Title
#'
#' @param x matrix of explanatory variables for log-relative risk
#' @param z matrix of explanatory variables for the nuisance parameter
#' @param y binary response variable
#' @param t binary risk factor
#' @param maxit maximum number of iteration for the estimation algorithm
#' @param start initial guess for the parameter. If null vector of zeros
#' @param param nuisance parameter. Possible values are "Richardson" or "Alternative"
#' @param method estimation method. Possible values are "Mle", "MeanBr" and "MedianBr"
#'
#' @return
#' @export
#'
#' @examples
brrr <- function(x,z=NULL,y,t,maxit=150, start=NULL,param = "Richardson", method= "Mle")
{

  if(is.null(z)) # this control was necessary for the the plot but it makes
    # redudant the controls in the estimating functions
  {
    z <- x
  }
  if(param == "Richardson")
  {
    est.function <- meth_richardson
  }
  else
  {
    est.function <- meth_alternative
  }

  mod0 = try(est.function(x,z,y,t,maxit=maxit, start=start, method= method))
  if(is.character(mod0))
  {
    mod0$convergence <- 0
  }
  if(mod0$convergence == 0)
  {
    mod0 = try(est.function(x,z,y,t,maxit=3*maxit, start=NULL, method= method, slowit = 1/2))
    if(is.character(mod0))
    {
      mod0$convergence <- 0
    }
  }

  if(method != "Mle")
  {
    if(mod0$convergence==1){
      mod1 <- try(est.function(x,z,y,t,maxit=4*maxit, start=mod0$Point.est, method= method,
                               slowit = 1/2))
      if(is.character(mod1)) mod1$convergence = 0
      if(mod1$convergence!=1)
      {
        mod1 <- try(est.function(x,z,y,t,maxit=4*maxit,start=mod0$Point.est,
                                 method= method,slowit = 1/4))
      }
    }
    else{
      mod1 <- try(est.function(x,z,y,t,maxit=maxit, method= method,slowit = 1/2))
      if(is.character(mod1)) mod1$convergence <- 0
      if(mod1$convergence!=1)
      {
        mod1 <- try(est.function(x,z,y,t,maxit=4*maxit, method= method, slowit = 1/4))
      }
    }
    mod0 <- mod1
  }
  # Quantities of interest
  pvalues <- 2*(1-stats::pnorm(abs(mod0$Point.est)/mod0$se))
  lower95 <- mod0$Point.est - mod0$se*stats::qnorm(0.975)
  upper95 <- mod0$Point.est + mod0$se*stats::qnorm(0.975)
  out <- data.frame(cbind( round(mod0$Point.est,4), round(mod0$se,4),
                           round(lower95 ,4),round(upper95 ,4), round(pvalues,4)))
  # names of the rows for out
  x <- as.matrix(x) # control for p = 1
  z <- as.matrix(z)
  rownames <- "gamma 0 "
  if(dim(x)[2]>1)
  {
    for(i in  1:(dim(x)[2]-1))
    {
      rownames <- c(rownames, paste("gamma",i,""))
    }
  }
  for(i in  0:(dim(z)[2]-1))
  {
    rownames <- c(rownames, paste("beta",i,""))
  }
  colnames(out) <- c("point est", "se", "lower 95%", "upper 95%", "p-value" )
  row.names(out) <- rownames
  output <- list(out = out, convergence = mod0$convergence, param = param,
                 method = method )
  return( structure(output, class = c("brrr")))
}
