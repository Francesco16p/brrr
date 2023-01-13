#' Improved and computationally stable relative risk estimation
#' for a binary exposure
#'
#' @param x matrix of explanatory variables for the nuisance parameter,
#' the first column is the intercept of the model
#' @param z matrix of explanatory variables for log-relative risk,
#' the first column is the intercept of the model
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
#' #####################################################
#' ## Relative risk regression with brrr            ####
#' #####################################################
#'
#'library(brrr)
#'library(HSAUR3)
#'
#' ##################################################
#' ###       Respiratory dataset                  ###
#' ### available on HSAUR3 R package              ###
#' ##################################################
#'
#'data <- HSAUR3::respiratory
#'
#'## Creating response variable, binary risk factor and model matrices
#'id <- seq.int(from = 5, to =  555, by  = 5)
#'data_matrix <- model.matrix(~., data[id, c(1,2,3,4,5)])
#'y <- data_matrix[,6]
#'t <- data_matrix[,3]
#'x <- z <- data_matrix[,-c(3,6)]
#'
#'## maximum likelihood fit using R package brrr, model: Richardson et al. (2017)
#'rr_mle <- brrr(x,z,y,t)
#'summary(rr_mle)
#'
#'## Mean bias-reduced fit using the R package brrr, model: Richardson et al. (2017)
#'rr_meanBR <- brrr(x,z,y,t, method = "MeanBr")
#'summary(rr_meanBR)
#'
#'## Median bias-reduced fit using the R function brrr, model: Richardson et al. (2017)
#'rr_medianBR <- brrr(x,z,y,t, method = "MedianBr")
#'summary(rr_medianBR)
#'
#'## maximum likelihood fit using R package brrr, model: Alternative
#'rr_mleA <- brrr(x,z,y,t, param = "Alternative")
#'summary(rr_mleA)
#'
#'## Mean bias-reduced fit using the R package brrr, model: Alternative
#'rr_meanBRA <- brrr(x,z,y,t, method = "MeanBr",param = "Alternative")
#'summary(rr_meanBRA)
#'
#'## Median bias-reduced fit using the R function brrr, model: Alternative
#'rr_medianBRA <- brrr(x,z,y,t, method = "MedianBr",param = "Alternative")
#'summary(rr_medianBRA)

brrr <- function(x, z=NULL, y, t, maxit=150, start=NULL, param = "Richardson", method= "Mle")
{

  # Change z and x so their are consistent with the notation in Kenne Pagui et al. (2023)
  x_swap <- z
  z_swap <- x

  x <- x_swap
  z <- z_swap

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

  #if(method != "Mle")
  #{
  #  if(mod0$convergence==1){
  #    mod1 <- try(est.function(x,z,y,t,maxit=4*maxit, start=mod0$Point.est, method= method,
  #                             slowit = 1/2))
  #    if(is.character(mod1)) mod1$convergence = 0
  #    if(mod1$convergence!=1)
  #    {
  #      mod1 <- try(est.function(x,z,y,t,maxit=4*maxit,start=mod0$Point.est,
  #                               method= method,slowit = 1/4))
  #    }
  #  }
  #  else{
  #    mod1 <- try(est.function(x,z,y,t,maxit=maxit, method= method,slowit = 1/2))
  #    if(is.character(mod1)) mod1$convergence <- 0
  #    if(mod1$convergence!=1)
  #   {
  #      mod1 <- try(est.function(x,z,y,t,maxit=4*maxit, method= method, slowit = 1/4))
  #    }
  #  }
  #  mod0 <- mod1
  #}
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
