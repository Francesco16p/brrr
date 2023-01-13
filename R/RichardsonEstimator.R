meth_richardson <- function(x,z=NULL,y,t,ep=1e-8,maxit=200, start=NULL, method= "Mle",slowit =1)
{
  # control z
  if(is.null(z)) z <- x

  # transform x,z in matrix
  x <- as.matrix(x)
  z <- as.matrix(z)
  #print(z)
  # re-ordering the data
  if(dim(x)[2]>1)
  {
    x <- rbind(x[t==0,], x[t==1,])
  }
  else
  {
    x <- matrix(c(x[t==0,], x[t==1,]), ncol = 1)
  }
  if(dim(z)[2]>1)
  {
    z <- rbind(z[t==0,], z[t==1,])
  }
  else
  {
    z <- matrix(c(z[t==0,], z[t==1,]), ncol = 1)
  }
  y <- c(y[t==0], y[t==1])
  t <- c(t[t==0], t[t==1])

  # getting indicies
  nobs <- nrow(x)
  nex <- nobs-sum(t) # number of exposed units
  nvarx <- ncol(x)
  nvarz <- ncol(z)
  if(is.null(start))
  {
    start<- rep(0,nvarx + nvarz)
  }

  # define x_delta and z_delta
  #x_delta <- array(NA, dim = c(n,n,p))
  #for(i in 1:nvarx) x_delta[,,i] <- diag(x[,i])
  #z_delta <- array(NA, dim = c(n,n,q))
  #for(i in 1:nvarz) z_delta[,,i] <- diag(z[,i])

  # defining functions for probabilities
  # function for Richardson's nusiance parameter
  prob <- function(eta1,eta2)
  {
    p0 <- rep(-1,length(eta1))
    num <- ( - (exp(eta1) + 1) * exp(eta2) + sqrt(exp(2 * eta2) * (exp(eta1) +
                                                                     1) ^ 2 + 4 * exp(eta1 + eta2) * (1 - exp(eta2))))
    den <- (2 * exp(eta1) * (1 - exp(eta2)))
    p0 <-  num/den
    p0[abs(eta2)< ep] <- 1 / (1 + exp(eta1[abs(eta2)< ep]))
    # 1e-8 it's an arbitrary value add a controll in next versions
    return(p0)
  }


  # derivatives of prob.standard
  # derivative of pi0 with respect to eta1
  p0_dp0_eta1 <- function(eta1, eta2)
  {
    dp0.eta1 <- ( -exp(eta2 - eta1) / (2 * (exp(eta2) - 1)) + exp(eta2) /
                    ((exp(eta2) - 1) * sqrt(4 * exp(eta2 + eta1) + (exp(eta1) -
                                                                      1)^2 *  exp(2 * eta2))) + (exp( - eta1) * (1 - exp(eta1)) *
                                                                                                   exp(2 * eta2)) / (2 * (exp(eta2) - 1)*sqrt(4 * exp(eta2 +
                                                                                                                                                        eta1) + (exp(eta1) - 1)^2 * exp(2*eta2))) )
    # extension for continuity
    dp0.eta1[abs(eta2)< ep] <- (- exp(eta1[abs(eta2)< ep]) /
                                  (exp(eta1[abs(eta2)< ep]) + 1) ^ 2)
    return(dp0.eta1)
  }
  # derivative of pi0 with respect to eta2
  # written in a different way compared to thesis
  p0_dp0_eta2 <- function(eta1,eta2)
  {
    dp0.eta2 <- ( - ((exp(eta1) + 1) * exp(eta2)) / (2 * exp(eta1) *
                   (exp(eta2) - 1) ^ 2) + exp(eta2) / ((exp(eta2) - 1) ^ 2 *
                    sqrt(4 * exp(eta2 + eta1) + (exp(eta1) - 1) ^ 2 * exp(2 * eta2))
                   ) + (exp( - eta1) * (exp(2 * eta1) + 1) * exp(2 * eta2)) /
                    (2 * (exp(eta2) - 1) ^ 2 * sqrt(4 * exp(eta2 + eta1) + (
                    exp(eta1) - 1) ^ 2 * exp(2 * eta2))))
    # extension for continuity
    dp0.eta2[abs(eta2)< ep] <- (exp(eta1[abs(eta2)< ep]) /
                                  (exp(eta1[abs(eta2)< ep]) + 1) ^ 3)
    return(dp0.eta2)
  }
  # Second derivative of pi0 with respect to eta1
  p0_d2p0_eta1 <- function(eta1,eta2)
  {
    d2p0.eta1 <- (exp(eta2 - eta1) / (2 * (exp(eta2) - 1)) + (exp(2 * eta2) *
                                                                (exp(eta2 + eta1) - exp(eta2) + 2) * (exp(eta2 + eta1) -
                                                                                                        exp(eta2) - 2 * exp(eta1))) / (2 * (exp(eta2) - 1) * (4 *
                                                                                                                                                                exp(eta2 + eta1)+(exp(eta1) - 1) ^ 2 * exp(2 * eta2)) ^ (3 / 2)
                                                                                                        ) + ( - (1 - exp(eta1)) * exp(2 * eta2 - eta1) - exp(2 * eta2)
                                                                                                        ) / (2 * (exp(eta2) - 1) * sqrt(4 * exp(eta2 + eta1) +
                                                                                                                                          (exp(eta1) - 1) ^ 2 * exp(2 * eta2))))
    # extension for continuity
    d2p0.eta1[abs(eta2)< ep] <-(((exp(eta1[abs(eta2)< ep]) - 1) *
                                   exp(eta1[abs(eta2)< ep])) / (exp(
                                     eta1[abs(eta2)< ep]) + 1) ^ 3)
    return(d2p0.eta1)
  }
  # Second derivative of pi0 with respect to eta2
  p0_d2p0_eta2 <- function(eta1,eta2)
  {
    d2p0.eta2 <- (((exp(eta1) + 1) * exp(2 * eta2-eta1)) / (exp(eta2) - 1) ^ 3
                  - ((exp(eta1) + 1) * exp(eta2 - eta1)) / (2 * (exp(eta2) - 1)
                                                            ^ 2) + ( - exp(2 * eta2) - exp(eta2)) / ((exp(eta2) - 1) ^ 3 *
                                                                                                       sqrt(4 * exp(eta2 + eta1) + (exp(eta1) - 1) ^ 2 * exp(2 *eta2))
                                                            ) + -((exp(2 * eta1) + 1) * exp(2 * eta2 - eta1)) /
                    ((exp(eta2) - 1) ^ 3 * sqrt(4 * exp(eta2 + eta1) + (exp(eta1) -
                                                                          1) ^ 2 * exp(2 * eta2))) -((exp(2 * eta1) + 1) * (4 * exp(eta2 +
                                                                                                                                      eta1) + 2 * (exp(eta1) - 1) ^ 2 * exp(2 * eta2)) * exp(2 *
                                                                                                                                                                                               eta2 - eta1)) / (4 * (exp(eta2) - 1) ^ 2 * (4 * exp(eta2 +
                                                                                                                                                                                                                                                     eta1) + (exp(eta1) - 1) ^ 2 * exp(2 * eta2)) ^ (3 / 2)) +
                    - (exp(eta2) * (4 * exp(eta2 + eta1) + 2 * (exp(eta1) - 1)
                                    ^ 2 *exp(2*eta2))) / (2 * (exp(eta2) - 1) ^ 2 * (4 * exp(eta2 +
                                                                                               eta1) + (exp(eta1) - 1) ^ 2 * exp(2 * eta2)) ^ (3 / 2)))
    # extension for continuity
    d2p0.eta2[abs(eta2)< ep] <- ( - ((exp(eta1[abs(eta2)< ep]) - 1) ^ 2 * exp(
      eta1[abs(eta2)< ep])) / (exp(
        eta1[abs(eta2)< ep]) + 1) ^ 5)
    return(d2p0.eta2)
  }
  # Second derivative of p0 respect to eta1 and eta2
  p0_d2p0_eta1_eta2 <- function(eta1,eta2)
  {
    d2p0.eta12 <- (exp(2 * eta2-eta1) / (2 * (exp(eta2) - 1) ^ 2) - exp(eta2 -
                                                                          eta1) / (2 * (exp(eta2) - 1)) + ((exp(eta1) - 1) * exp(3 *
                                                                                                                                   eta2 - eta1) - 2 * exp(eta2)) / (2 * (exp(eta2) - 1) ^ 2 *
                                                                                                                                                                      sqrt(4 * exp(eta2 + eta1) + (exp(eta1) - 1) ^ 2 * exp(2 *
                                                                                                                                                                                                                              eta2))) + (1 - exp(eta1)) * exp(2 * eta2 - eta1) /
                     ((exp(eta2) - 1) * sqrt(4 * exp(eta2 + eta1) + (exp(eta1) - 1)
                                             ^ 2*exp(2*eta2))) - (exp(eta2) * (4 * exp(eta2 + eta1) + 2 *
                                                                                 (exp(eta1) - 1) ^ 2 * exp(2 * eta2))) / (2 * (exp(eta2) - 1) *
                                                                                                                            (4 * exp(eta2 + eta1) + (exp(eta1) - 1) ^ 2 * exp(2 * eta2))
                                                                                                                          ^ (3 / 2)) - ((1 - exp(eta1)) * (4 * exp(eta2 + eta1) + 2 *
                                                                                                                                                             (exp(eta1) - 1) ^ 2 * exp(2 * eta2)) * exp(2 * eta2 - eta1)) /
                     (4 * (exp(eta2) - 1) * (4 * exp(eta2 + eta1) + (exp(eta1) - 1)
                                             ^ 2 *exp(2 * eta2)) ^ (3 / 2)))
    # extension for continuity
    d2p0.eta12[abs(eta2)< ep] <- ( - (exp(eta1[abs(eta2)< ep]) * (2 *
                                                                    exp(eta1[abs(eta2)< ep]) - 1)) /
                                     (exp(eta1[abs(eta2)< ep]) + 1) ^ 4)
    return(d2p0.eta12)
  }

  # Key quantities
  key_quantities <- function(pars)
  {
    # Basic quantities
    gammas <- pars[1:nvarx]
    betas <- pars[-c(1:nvarx)]
    eta1 <- drop(x %*% gammas)
    eta2 <- drop(z %*% betas)
    pi <- prob(eta1,eta2) * exp(eta1 * t)
    d.1.eta1 <- p0_dp0_eta1(eta1,eta2) * exp(eta1 * t) + pi * t
    d.1.eta2 <- p0_dp0_eta2(eta1,eta2) * exp(eta1 * t)
    d.2.eta1 <- (p0_d2p0_eta1(eta1,eta2) * exp(eta1 * t) + 2 *
                   p0_dp0_eta1(eta1,eta2) * exp(eta1 * t) * t + pi * t)
    d.2.eta2 <- p0_d2p0_eta2(eta1,eta2) * exp(eta1 * t)
    d.2.eta12 <- (p0_d2p0_eta1_eta2(eta1,eta2) +
                    p0_dp0_eta2(eta1,eta2) * t) * exp(eta1 * t)

    # Diagonal matrices

    # d.1.eta1 <- diag(d.1.eta1)
    #d.1.eta2 <- diag(d.1.eta2)
    #d.2.eta1 <- diag(d.2.eta1)
    #d.2.eta2 <- diag(d.2.eta2)
    #d.2.eta12 <- diag(d.2.eta12)

    # Quanties related to the variance

    v <- pi * (1 - pi) # variance
    d1v <- 1 - 2 * pi  # variance's derivative w.r.t pi
    #V <- diag(v)
    #d1v <- diag(d1v)

    out <- list(betas = betas, gammas = gammas, eta1 = eta1, eta2 = eta2,
                pi = pi, d.1.eta1 = d.1.eta1, d.1.eta2 = d.1.eta2,
                d.2.eta1 = d.2.eta1, d.2.eta2 = d.2.eta2,  d.2.eta12 = d.2.eta12,
                v = v, d1v = d1v)
    return(out)
  }
  # Score-function
  score <- function(pars,fit=NULL)
  {
    if(is.null(fit))
    {
      fit <- key_quantities(pars)
    }

    out <-with( fit,
                # score for each observation
                {
                  score_gamma <- ( x * d.1.eta1 * (y-pi) / v)
                  score_beta <-  ( z * d.1.eta2 * (y-pi) / v)
                  score_tot <- colSums(cbind(score_gamma, score_beta))
                  score_tot
                }
    )
    out
  }

  # Fisher expected information
  information <- function(pars,inverse = FALSE,fit = NULL)
  {
    if(is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    info<- with(fit,{
      # defining the blocks of the expected information
      igg <- crossprod( x , (d.1.eta1 ^ 2 * x / v))
      igb <- crossprod( x , (d.1.eta1 * d.1.eta2 * z / v))
      ibb <- crossprod( z , (d.1.eta2 ^ 2 * z / v))
      info <- as.matrix(rbind( cbind(igg, igb),
                               cbind(t(igb), ibb)))
      info
    })
    if(!inverse)
    {
      return(info)
    }
    else
    {
      return(chol2inv(chol(info)))
    }
  }

  # nvarx + nvarz gamma matrix
  P_Q_gamma <- function(pars, fit = NULL)
  {
    if(is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    PQg <- with(fit,{
      # parts costant in the matrix
      PQg <- array(data = NA, dim= c(nvarx + nvarz, nvarx + nvarz, nvarx))
      for(t in 1:nvarx)
      {
        xt =x[, t]
        L.gg.t <- crossprod( x ,( d.2.eta1 * d.1.eta1 / v ) * xt * x)
        L.gb.t <- crossprod( x ,( d.2.eta12 * d.1.eta1 / v ) * xt * z)
        L.bb.t <- crossprod( z ,( d.2.eta2 * d.1.eta1 / v  ) * xt * z)

        PQg[,,t] <- as.matrix( rbind( cbind(L.gg.t,L.gb.t),
                                      cbind(t(L.gb.t),L.bb.t ) ))


      }
      PQg
    })
    return(PQg)
  }

  # P+ Q beta matrix
  P_Q_beta <- function(pars, fit = NULL)
  {
    if(is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    PQb <- with(fit,{
      # parts costant in the matrix
      PQb <- array(data = NA, dim= c(nvarx + nvarz,nvarx + nvarz, nvarz))

      for(s in 1:nvarz)
      {
        zs <- z[,s]
        L.gg.s <- crossprod( x ,( d.2.eta1 * d.1.eta2 / v )*zs * x)
        L.gb.s <- crossprod( x ,( d.2.eta12 * d.1.eta2 / v ) * zs * z)
        L.bb.s <- crossprod( z ,( d.2.eta2 * d.1.eta2 / v) * zs * z)
        PQb[,,s] <- as.matrix( rbind( cbind(L.gg.s,L.gb.s),
                                      cbind(t(L.gb.s),L.bb.s)))
      }
      PQb
    })
    return(PQb)
  }

  # P/3+Q/2 gamma matrix
  P_Q_3_2_gamma <- function(pars, fit = NULL)
  {
    if(is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    PQ32g <- with(fit,{
      # parts costant in the matrix
      PQ32g <- array(data = NA, dim= c(nvarx + nvarz,nvarx + nvarz, nvarx))

      for(t in 1:nvarx)
      {
        xt <- x[,t]
        L.gg.t <- crossprod( x ,( 1 / 3 * d1v * d.1.eta1 ^ 3 / v ^ 2 +
                                    1 / 2 * (d.2.eta1 * d.1.eta1 / v - d1v *
                                               d.1.eta1 ^ 3 / v ^ 2)) * xt * x)
        L.gb.t <- crossprod( x ,( 1 / 3 * d1v * d.1.eta1 ^ 2 * d.1.eta2 / v ^ 2 +
                                    1 / 2 * ( d.2.eta12 * d.1.eta1 / v - d1v *
                                                d.1.eta1 ^ 2 * d.1.eta2 / v ^ 2)) * xt * z)
        L.bb.t <- crossprod( z ,( 1 / 3 * d1v * d.1.eta1 * d.1.eta2 ^ 2 / v ^ 2 +
                                    1 / 2 * (d.2.eta2 * d.1.eta1 / v - d1v *
                                               d.1.eta1 * d.1.eta2 ^ 2 / v ^ 2)) * xt * z)

        PQ32g[,,t] <- as.matrix( rbind( cbind(L.gg.t,L.gb.t),
                                        cbind(t(L.gb.t),L.bb.t ) ))
      }
      PQ32g
    })
    return(PQ32g)
  }

  # P/3+ Q/2 beta matrix
  P_Q_3_2_beta <- function(pars, fit = NULL)
  {
    if(is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    PQ32b <- with(fit,{
      # parts constant in the matrix
      PQ32b <- array(data = NA, dim= c(nvarx + nvarz,nvarx + nvarz,nvarz))

      for(s in 1:nvarz)
      {
        zs <- z[,s]
        L.gg.s <- crossprod( x ,(1 / 3 * d1v * d.1.eta1 ^ 2 * d.1.eta2 / v ^ 2 +
                                   1 / 2 * (d.2.eta1 * d.1.eta2 / v - d1v *
                                              d.1.eta1 ^ 2 * d.1.eta2 / v ^ 2)) * zs * x)
        L.gb.s <- crossprod( x ,(1 / 3 * d1v * d.1.eta1 * d.1.eta2 ^ 2 / v ^ 2 +
                                   1 / 2 * (d.2.eta12 * d.1.eta2 / v - d1v *
                                              d.1.eta1 * d.1.eta2 ^ 2 / v ^ 2)) * zs * z)
        L.bb.s <- crossprod( z ,(1 / 3 * d1v * d.1.eta2 ^ 3 / v ^ 2 +
                                   1 / 2 * (d.2.eta2 * d.1.eta2 / v - d1v *
                                              d.1.eta2 ^ 3 / v ^ 2)) * zs * z)

        PQ32b[,,s] <- as.matrix( rbind( cbind(L.gg.s,L.gb.s),
                                        cbind(t(L.gb.s),L.bb.s ) ))
      }
      PQ32b
    })
    return(PQ32b)

  }
  # Matrix F1
  F1 <- function(pars, fit = NULL)
  {
    if(is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    f1 <- rep(NA,nvarx + nvarz)
    p.q.gamma <- P_Q_gamma(fit=fit)
    p.q.beta <- P_Q_beta(fit=fit)
    info.inv <- information(inverse = TRUE, fit = fit)
    for(k in 1:(nvarx+nvarz))
    {
      if( k<=nvarx)
      {
        f1[k] <- sum(diag(info.inv %*% p.q.gamma[,, k])) # maybe change %*%
      }
      else
      {
        f1[k] <- sum(diag(info.inv %*% p.q.beta[,, k-nvarx] ))
      }

    }
    return(f1)
  }

  # h function
  h_function <- function(pars, fit = NULL)
  {
    if(is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    h <- with(fit,{
      h <- array(data = NA, dim= c(nvarx + nvarz, nvarx + nvarz, nvarx + nvarz))
      info.inv <- information(inverse= TRUE, fit=fit)
      for( k in 1:(nvarx + nvarz) )
      {
        h[, , k] <-  tcrossprod( info.inv[, k], info.inv[, k]) / info.inv[k, k]
      }
      h
    })
    return(h)
  }

  # Creating F2 a function that calculate all the 2*(p-1) F2 vectors


  # Creating F2 a function that calculate all the 2*(p-1) F2 vectors
  F2_tilde <- function(pars, fit = NULL)
  {
    if(is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    f2 <- matrix(NA,ncol = nvarx + nvarz, nrow = nvarx + nvarz )
    f2.tilde <- rep(NA, nvarx + nvarz)
    p.q.3.2.gamma <- P_Q_3_2_gamma(fit = fit)
    p.q.3.2.beta <- P_Q_3_2_beta(fit = fit)
    h<- h_function(fit=fit)
    inv.info <- information(inverse = TRUE, fit = fit)
    for(k in 1:(nvarx + nvarz))
    {
      if(k <= nvarx){
        for(j in 1:(nvarx + nvarz)){
          f2[j, k] <- sum(diag(h[, , j] %*% p.q.3.2.gamma[, , k]))
        }
      }
      else
      {
        for(j in 1:(nvarx + nvarz))
        {
          f2[j,k] <- sum(diag(h[,,j]%*%p.q.3.2.beta[,, k - nvarx ]))
        }
      }

    }
    for(l in 1:(nvarx + nvarz))
    {
      f2.tilde[l] <- inv.info[l, ]%*%f2[l,]
    }
    return(f2.tilde)
  }
  ## median  adjustment term


  #Estimation
  max.step.factor<-12
  epsilon <- 1e-5 #  qua hai cambiato
  par <- start ## initial value
  # Controll for mle, mean bias reduction and median
  # bias reduction
  controlA <- 0
  controlB <- 0
  # Mean bias reduction
  if(method=="MeanBr")
  {
    controlA <- 1
  }
  # Median bias reduction
  if(method=="MedianBr")
  {
    controlA <- 1
    controlB <- 1
  }
  fit <- key_quantities(par)
  adjusted.grad <- (score(fit=fit) + 1/2 * F1(fit=fit) * controlA -
                      information(fit=fit) %*%  F2_tilde(fit=fit) * controlB)
  step.par <- information(fit=fit, inverse = TRUE) %*% adjusted.grad
  for (iter in seq.int(maxit))
  {
    step.factor <- 0
    testhalf <- TRUE
    while (testhalf & step.factor < max.step.factor)
    {
      step.par.previous <- step.par
      par <- par + slowit * 2 ^ ( - step.factor) * step.par
      fit <- key_quantities(par)
      adjusted.grad <- (score(fit=fit) + 1 / 2 * F1(fit=fit) * controlA -
                          information(fit=fit) %*% F2_tilde(fit=fit) * controlB)
      step.par <- information(fit=fit, inverse = TRUE) %*% adjusted.grad
      if (step.factor == 0 & iter == 1)
      {
        testhalf <- TRUE
      }
      else
      {
        testhalf <- sum(abs(step.par),na.rm=TRUE) > sum(abs(step.par.previous),na.rm = TRUE)
      }
      step.factor <- step.factor + 1
    }


    if ( (sum(abs(adjusted.grad),na.rm = TRUE) < epsilon) | ( any(abs(par)>15) ) ) {
      break # was (sum(abs(step.par),na.rm = TRUE) < epsilon)
    }
  }


  if(iter >= maxit| (any(abs(par) > 15) > 0))
  {
    convergence <- 0
    warning("optimization failed to converge")
  }
  else
  {
    convergence <- 1
  }
  Point.est <- par
  var.est <- information(Point.est, inverse = TRUE)
  se <- sqrt(diag(var.est))
  return(list(Point.est=Point.est, se=se,convergence=convergence))

}
