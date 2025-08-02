#' PQL function for binary data
#'
#' @keywords internal
binomial_PQL <- function(Y,X_sig1=NULL, Beta, Z, Alpha){
  if(is.null(X_sig1)){
    exp_value <- exp(Beta+ as.matrix(rowSums(mapply(function(x,y){x%*%y}, Z,Alpha)), ncol=1))
  }else{
    exp_value <- exp(X_sig1%*%Beta+ as.matrix(rowSums(mapply(function(x,y){x%*%y}, Z,Alpha)), ncol=1))
  }
  mu <- exp_value/(1+exp_value)
  V <- drop(exp_value/(1+exp_value)^2)
  inv_V_vector <- 1/V
  y_star <- inv_V_vector*(Y-mu)+log(exp_value)
  return(list(mu=mu, V=V, inv_V_vector=inv_V_vector, y_star=y_star))
}

#' PQL function for count data
#'
#' @keywords internal
poisson_PQL <- function(Y,X_sig1=NULL, Beta, Z, Alpha, offset){
  if(is.null(offset)){
    offset = 1
  }

  if(is.null(X_sig1)){
    mu <- exp(Beta+ as.matrix(rowSums(mapply(function(x,y){x%*%y}, Z,Alpha)), ncol=1))*offset
  }else{
    mu <- exp(X_sig1%*%Beta+ as.matrix(rowSums(mapply(function(x,y){x%*%y}, Z,Alpha)), ncol=1))*offset
  }

  V <- drop(mu)
  inv_V_vector <- 1/V
  y_star <- inv_V_vector*(Y-mu)+log(mu)-log(offset)
  return(list(mu=mu, V=V, inv_V_vector=inv_V_vector, y_star=y_star))
}

#' PQL function
#'
#' @keywords internal
PQL <- function(Y, Z, kinship, X=NULL, Xc=NULL, Xs=NULL, indices_X=NULL, indices_Xc=NULL, family, offset=NULL, postprob=NULL){
  # kinship: covariance matrix for random effects
  # X: SNPs from screening
  # Xc: SNPs from last ite
  # Xs: other covariates

  # family & default link function
  # binomial	(link = "logit")
  # poisson	(link = "log")

  if((is.null(X) + is.null(Xc) + is.null(Xs)) < 3) {
    if(!is.null(X) | !is.null(Xc)){
      X_sig <- cbind(X,Xc)
      indices <- c(indices_X, indices_Xc)

      if(ncol(X_sig)>1){
        if(Matrix::rankMatrix(X_sig)[1] < ncol(X_sig)){
          dropped_cols <- caret::findLinearCombos(X_sig)$remove
          X_sig <- X_sig[,-dropped_cols]
          indices <- indices[-dropped_cols]
        }

      }

      if(ncol(X_sig)>1){
        total.p <- ncol(X_sig)
        redundancy.snp.indices <- c()
        for(snpi in 1:(total.p-1)){
          for(snpj in (snpi+1):total.p){
            if(!(snpi %in% redundancy.snp.indices)){
              if(stats::cor.test(X_sig[,snpi],X_sig[,snpj],method = "spearman",exact=FALSE,alternative = "greater")$p.value < 0.05){
                if(postprob[indices[snpi]] < postprob[indices[snpj]]){
                  redundancy.snp.indices <- c(redundancy.snp.indices,snpi)
                }else{
                  redundancy.snp.indices <- c(redundancy.snp.indices,snpj)
                }
              }
            }
          }
        }
        X_sig <- X_sig[,!((1:total.p) %in% redundancy.snp.indices)]
        indices <- indices[!((1:total.p) %in% redundancy.snp.indices)]

      }

      if(is.null(Xs)){
        X_sig1 <- cbind(1,X_sig)
      }else{
        X_sig1 <- cbind(1,Xs,X_sig)
      }


    }else{
      X_sig1 <- cbind(1,Xs)
    }

    glmfit <- stats::glm(Y~X_sig1[,-1], family = family)
    rm(X)
    Beta <- glmfit$coefficients
    n <- length(Y)
    n_rf <- length(kinship)
    Alpha <- list()
    Kappa <- list()
    Kappa_temp <- list()
    for(i_rf in 1:n_rf){
      Alpha[[i_rf]] <- matrix(rep(0, nrow(kinship[[i_rf]])), ncol = 1)
      Kappa[[i_rf]] <- 0
      Kappa_temp[[i_rf]] <- 0
    }
    Beta_temp <- matrix(0, nrow = ncol(X_sig1), ncol = 1)


    while(sum(mapply(function(x,y){abs(x-y)>0.0001},Kappa, Kappa_temp))>=1 | max(abs(Beta-Beta_temp))>0.0001){

      Kappa_temp <- Kappa
      Beta_temp <- Beta

      if(family=="binomial"){
        par_est <- binomial_PQL(Y=Y,X_sig1=X_sig1, Beta=Beta, Z=Z, Alpha=Alpha)
        mu <- par_est$mu
        V <- par_est$V
        inv_V_vector <- par_est$inv_V_vector
        y_star <- par_est$y_star

        V_05 <- sqrt(V)
        VkV <- V_05*matrix(mapply(function(y,z){z%*%y%*%t(z)}, kinship, Z), ncol=n, nrow=n)%*%diag(V_05)

        eign <- eigen(VkV, symmetric = TRUE)
        P <- eign$vectors
        D <- eign$values
        D[D<10^(-10)] = 0

        re <- y_star-X_sig1%*%Beta
        reP <- matrix(re*V_05,nrow=1)%*%P    # re^T %*% V^0.5 %*%P
        XP <- t(V_05*X_sig1) %*% P  # X^T %*% V^0.5 %*%P
        logl <- function(x){
          # log posterior of kappa
          XHX <- XP %*% ((x*D+1)*t(XP))
          l <- (-0.5*sum(log(inv_V_vector))-0.5*sum(log(x*D+1)) #log(det(H))
                -0.5*determinant(XHX, logarithm = T)$modulus[1]  # t(X) %*% inv H %*% X
                -0.5*sum(reP^2/(x*D+1))  #t(re) %*% inv H %*% re = t(re.tilde) %*% inv D %*% re.tilde
          )
          return(l)
        }

        Kappa <- stats::optim(par=0.1, fn=logl, lower =0 , upper = 5, method = "L-BFGS-B", control = list(fnscale=-1))$par
        Kappa <- as.list(Kappa)

      }

      if(family=="poisson"){
        par_est <- poisson_PQL(Y=Y,X_sig1=X_sig1, Beta=Beta, Z=Z, Alpha=Alpha, offset = offset)
        mu <- par_est$mu
        V <- par_est$V
        inv_V_vector <- par_est$inv_V_vector
        y_star <- par_est$y_star

        logl <- function(x){
          # log posterior of kappa
          H <- matrix(rowSums(mapply(function(x,y,z){x*z%*%y%*%t(z)}, x, kinship, Z)), ncol=n, nrow=n)
          diag(H) <- diag(H)+inv_V_vector
          eign <- eigen(H, symmetric = TRUE)
          P <- eign$vectors
          D <- eign$values
          D_inv <- 1/D
          D_inv[!is.finite(D_inv)] <- 0
          X.significant_1.tilde <- t(P)%*%X_sig1  # t(P) %*% (1 X)

          re.tilde <- t(P)%*%y_star-X.significant_1.tilde%*%Beta

          l <- (+0.5*sum(log(D_inv))   #log(det(H))
                -0.5*determinant(t(X.significant_1.tilde)%*%(D_inv*X.significant_1.tilde), logarithm=TRUE)$modulus[1] #t(X) %*% inv_H %*% X = t(X.tilde) %*% D_inv %*% X.tilde
                -0.5*sum(re.tilde^2*D_inv)) # (y-Xbeta)invH(y-Xbeta)
          return(l)
        }

        ts <- stats::optim(par=rep(0.1,n_rf), fn=logl, lower = rep(0,n_rf), upper = rep(50,n_rf), method = "L-BFGS-B", control = list(fnscale=-1))$par
        Kappa <- as.list(ts)

      }

      H <- matrix(rowSums(mapply(function(x,y,z){x*z%*%y%*%t(z)}, Kappa, kinship, Z)), ncol=n, nrow=n)
      diag(H) <- diag(H)+inv_V_vector
      eign <- eigen(H, symmetric = TRUE)
      D <- eign$values
      P <- eign$vectors
      D_inv <- 1/D
      D_inv[!is.finite(D_inv)] <- 0

      X.significant_1.tilde <- t(P)%*%X_sig1
      y.tilde <- t(P)%*%y_star


      #update beta and alpha
      xdx <- t(X.significant_1.tilde)%*%(D_inv*X.significant_1.tilde)
      eign <- eigen(xdx, symmetric = TRUE)
      D_xdx <- eign$values
      P_xdx <- eign$vectors
      D_inv_xdx <- 1/D_xdx
      D_inv_xdx[!is.finite(D_inv_xdx)] <- 0
      inv_xdx <- P_xdx%*%(D_inv_xdx*t(P_xdx))
      Beta <- inv_xdx%*%t(X.significant_1.tilde)%*%matrix(D_inv*y.tilde,ncol = 1)
      re.tilde <- y.tilde-X.significant_1.tilde%*%Beta   #t(P)%*%(y_star-XBeta)
      Alpha <- mapply(function(x,y,z){x*y%*%t(z)%*%P%*%matrix(D_inv*re.tilde,ncol=1)}, Kappa, kinship, Z, SIMPLIFY = FALSE)    #Kappa*kinship%*%t(Z)%*%inv_H%*%(y_star-XBeta)

    }


    if(family=="binomial"){
      par_est <- binomial_PQL(Y=Y,X_sig1=X_sig1, Beta=Beta, Z=Z, Alpha=Alpha)
    }
    if(family=="poisson"){
      par_est <- poisson_PQL(Y=Y,X_sig1=X_sig1, Beta=Beta, Z=Z, Alpha=Alpha, offset = offset)
    }


    mu <- par_est$mu
    V <- par_est$V
    inv_V_vector <- par_est$inv_V_vector
    y_star <- par_est$y_star

    return(list(y_star=y_star, kappa=Kappa, H=H, P=P, D_inv=D_inv, beta=Beta, inv_v=inv_V_vector, X_sig1=X_sig1))

  }else{

    glmfit <- stats::glm(Y~1, family = family)
    Beta <- glmfit$coefficients
    n <- length(Y)
    n_rf <- length(kinship)
    Alpha <- list()
    Kappa <- list()
    Kappa_temp <- list()
    for(i_rf in 1:n_rf){
      Alpha[[i_rf]] <- matrix(rep(0, nrow(kinship[[i_rf]])), ncol = 1)
      Kappa[[i_rf]] <- 0
      Kappa_temp[[i_rf]] <- 0
    }
    Beta_temp <- 1

    while(sum(mapply(function(x,y){abs(x-y)>0.0001},Kappa, Kappa_temp))>=1 | abs(Beta-Beta_temp)>0.0001){

      Kappa_temp <- Kappa
      Beta_temp <- Beta

      if(family=="binomial"){
        par_est <- binomial_PQL(Y=Y,X_sig1=NULL, Beta=Beta, Z=Z, Alpha=Alpha)
        mu <- par_est$mu
        V <- par_est$V
        inv_V_vector <- par_est$inv_V_vector
        y_star <- par_est$y_star

        V_05 <- sqrt(V)

        VkV <- V_05*matrix(mapply(function(y,z){z%*%y%*%t(z)}, kinship, Z), ncol=n, nrow=n)%*%diag(V_05)

        eign <- eigen(VkV, symmetric = TRUE)
        P <- eign$vectors
        D <- eign$values
        D[D<10^(-10)] = 0
        re <- y_star-Beta
        V05P <- matrix(V_05,nrow=1)%*%P # 1^T %*% (V^0.5) %*% P
        reP <- matrix(re*V_05,nrow=1)%*%P    # re^T %*% V^0.5 %*%P


        logl <- function(x){
          # log posterior of kappa
          l <- (-0.5*sum(log(inv_V_vector))-0.5*sum(log(x*D+1)) #log(det(H))
                -0.5*log(sum(V05P^2/(x*D+1))) # t(1) %*% inv H %*% 1 = t(1) %*% P %*% inv D %*% t(P) %*% 1
                -0.5*sum(reP^2/(x*D+1))  #t(re) %*% inv H %*% re = t(re.tilde) %*% inv D %*% re.tilde
          )
          return(l)
        }

        Kappa <- stats::optim(par=0.1, fn=logl, lower =0 , upper = 30, method = "L-BFGS-B", control = list(fnscale=-1))$par
        Kappa <- as.list(Kappa)
      }

      if(family=="poisson"){
        par_est <- poisson_PQL(Y=Y,X_sig1=NULL, Beta=Beta, Z=Z, Alpha=Alpha, offset = offset)
        mu <- par_est$mu
        V <- par_est$V
        inv_V_vector <- par_est$inv_V_vector
        y_star <- par_est$y_star

        logl <- function(x){
          # log posterior of kappa
          H <- matrix(rowSums(mapply(function(x,y,z){x*z%*%y%*%t(z)}, x, kinship, Z)), ncol=n, nrow=n)
          diag(H) <- diag(H)+inv_V_vector
          eign <- eigen(H, symmetric = TRUE)
          P <- eign$vectors
          D <- eign$values
          D_inv <- 1/D
          D_inv[!is.finite(D_inv)] <- 0
          one.tilde <- colSums(P)    # t(P) %*% 1
          #inv_H <- P%*%(D_inv*t(P))
          re.tilde <- t(P)%*%(y_star-Beta)

          l <- (+0.5*sum(log(D_inv)) #log(det(H))
                -0.5*log(sum(one.tilde^2*D_inv)) # t(1) %*% inv H %*% 1 = t(1) %*% P %*% inv D %*% t(P) %*% 1
                -0.5*sum(re.tilde^2*D_inv))  #t(re) %*% inv H %*% re = t(re.tilde) %*% inv D %*% tilde
          return(l)
        }

        ts <- stats::optim(par=rep(0.1,n_rf), fn=logl, lower = rep(0,n_rf), upper = rep(30,n_rf), method = "L-BFGS-B", control = list(fnscale=-1))$par
        Kappa <- as.list(ts)
      }

      H <- matrix(rowSums(mapply(function(x,y,z){x*z%*%y%*%t(z)}, Kappa, kinship, Z)), ncol=n, nrow=n)
      diag(H) <- diag(H)+inv_V_vector
      eign <- eigen(H, symmetric = TRUE)
      D <- eign$values
      P <- eign$vectors
      D_inv <- 1/D
      D_inv[!is.finite(D_inv)] <- 0


      X.tilde <- colSums(P)
      y.tilde <- t(P)%*%y_star
      #update beta and alpha
      xdx <- sum(X.tilde^2*D_inv)
      inv_xdx <- 1/xdx
      Beta <- inv_xdx*sum(X.tilde*D_inv*y.tilde)

      re.tilde <- y.tilde-X.tilde*Beta   #t(P)%*%(y_star-XBeta)
      Alpha <- mapply(function(x,y,z){x*y%*%t(z)%*%P%*%matrix(D_inv*re.tilde,ncol=1)}, Kappa, kinship, Z, SIMPLIFY = FALSE)
      #Kappa*kinship%*%t(Z)%*%inv_H%*%(y_star-XBeta)

    }
    if(family=="binomial"){
      par_est <- binomial_PQL(Y=Y,X_sig1=NULL, Beta=Beta, Z=Z, Alpha=Alpha)
    }

    if(family=="poisson"){
      par_est <- poisson_PQL(Y=Y,X_sig1=NULL, Beta=Beta, Z=Z, Alpha=Alpha, offset = offset)
    }

    mu <- par_est$mu
    V <- par_est$V
    inv_V_vector <- par_est$inv_V_vector
    y_star <- par_est$y_star

    return(list(y_star=y_star, kappa=Kappa, H=H, P=P, D_inv=D_inv, beta=Beta, inv_v=inv_V_vector))
  }

}

#' likelihood function for null
#'
#' @keywords internal
log_marginal_likelihood_null <- function(y.tilde,D_inv){
  # This function computes the log marginal likelihood for the case
  # when there is no regressor in the model. y.tilde ~ N(0,D)
  n <- length(y.tilde)
  return(-0.5*n*log(2*pi)+0.5*sum(log(D_inv))-0.5*sum(D_inv*y.tilde^2))
}

#' likelihood function
#'
#' @keywords internal
log.marginal.likelihood <- function(k, x.tilde_m, y.tilde, D_inv, ydinvy, dinvy, g)
{

  n <- length(y.tilde)     # sample size
  #k <- ncol(x.tilde_m)
  return( -0.5*n*log(2*pi) - 0.5*k*log(1+g) + 0.5*sum(log(D_inv))
          -0.5*ydinvy + 0.5*g/(g+1)*t(dinvy)%*% x.tilde_m %*% solve(t(x.tilde_m)%*%(D_inv*x.tilde_m)) %*% t(x.tilde_m) %*% dinvy
  )
}

#' GINAX function
#'
#' @keywords internal
GINAX_terminal <- function(Y, kinship, Z, SNPs, family, offset=NULL,
                           FDR.threshold, maxiterations, runs_til_stop){
  ite <- 1
  indices_sig_list <- list()
  Xc <- NULL
  X = NULL
  Xs=NULL
  pi0_1 <- NULL
  indices_Xc <- NULL
  tmp = NULL
  indices_previous = NULL
  postprob = NULL

  # screen_list = list()
  # select_list = list()

  while((ite <= 10) &
        ((is.null(tmp)) |
         ( (!is.character(tmp$modelselection))
           & !((length(indices_Xc) == length(indices_previous))
               & sum(!(indices_Xc %in% indices_previous))==0))) ){
    # stop criterion:
    # num of ite larger than 10
    # screening does not select any SNP
    # SNPs selected in this ite are the same as those in the last ite.

    indices_previous <- indices_Xc

    n = length(Y)
    PQL_est <- PQL(Y=Y, Z=Z, kinship=kinship, X=NULL, Xc=Xc, Xs=Xs, indices_X=NULL, indices_Xc=indices_Xc, family=family, postprob=postprob, offset=NULL)

    y_star <- PQL_est$y_star
    Kappa <- PQL_est$kappa
    Beta <- PQL_est$beta
    inv_V_vector <- PQL_est$inv_v
    H <- PQL_est$H
    P <- PQL_est$P
    D_inv <- PQL_est$D_inv

    if(is.null(Xc) & is.null(Xs)){

      A <- H-sum(P%*%(D_inv*t(P)))
      eign <- eigen(A, symmetric = TRUE)
      P <- eign$vectors
      y.tilde = t(P) %*% (y_star-Beta)
      D <- eign$values
      D_inv <- D^(-1)
      #D_inv[!is.finite(D_inv)] <- 0
      D_inv[n] <- 0
    }else{
      X_1 <-  PQL_est$X_sig1
      A <- H-X_1%*%solve(t(X_1)%*%P%*%(D_inv*t(P))%*%X_1)%*%t(X_1)
      eign <- eigen(A, symmetric = TRUE)
      P <- eign$vectors
      y.tilde = t(P) %*% (y_star-X_1%*%Beta)
      D <- eign$values
      D_inv <- D^(-1)
      D_inv[(n-ncol(X_1)+1):n] <- 0
    }

    X.tilde = t(P) %*% SNPs
    # Estimate beta with beta.hat and compute var(beta.hat)
    xj.t.xj = apply(X.tilde*(D_inv*X.tilde),2,sum)
    xj.t.y = t(X.tilde) %*% (D_inv*y.tilde)

    beta.hat = xj.t.y / xj.t.xj
    var.beta.hat = 1 / xj.t.xj
    t.statistic <- beta.hat / sqrt(var.beta.hat)
    pvalues <- 2*stats::pnorm(abs(t.statistic), mean=0, sd=1, lower.tail=FALSE)
    P3D_return_dat <- cbind(beta.hat,var.beta.hat,pvalues)
    P3D_return_dat <- as.data.frame(P3D_return_dat)
    colnames(P3D_return_dat) <- c("Beta_Hat","Var_Beta_Hat","P_Values")
    #indices_Xc=indices_previous
    n= length(Y)
    #estimate pi0
    log.marginal.likelihood_par = function(param)
    {
      pi0 = exp(param[1]) / (1+exp(param[1]))
      g = exp(param[2])
      return(sum(log(pi0*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat)) +
                       (1-pi0)*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat*(g+1)))
      )))

    }

    result <- stats::optim(c(2,-2), fn=log.marginal.likelihood_par, lower = c(1,-20), method = "L-BFGS-B", hessian=TRUE, control = list(fnscale=-1))

    pi0.hat = exp(result$par[1]) / (1+exp(result$par[1]))
    g.hat = exp(result$par[2])

    # Compute posterior probability of beta_j different than 0:
    numerator <- (1-pi0.hat)*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat*(g.hat+1)))
    denominator <- pi0.hat*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat)) +
      (1-pi0.hat)*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat*(g.hat+1)))

    postprob <- numerator / denominator

    if(!is.null(indices_previous)){
      postprob[indices_previous] = 0
    }

    order.postprob <- order(postprob, decreasing=TRUE)
    postprob.ordered <- postprob[order.postprob]


    FDR.Bayes <- cumsum(postprob.ordered) / 1:ncol(SNPs)
    if(sum(FDR.Bayes > FDR.threshold) == 0){
      P3D_return_dat <- cbind(P3D_return_dat,postprob,FALSE)
      P3D_return_dat <- as.data.frame(P3D_return_dat)
      colnames(P3D_return_dat) <- c("Beta_Hat","Var_Beta_Hat","P_Values","PostProb","Significant")
    }else{
      P3D_return_dat <- cbind(P3D_return_dat,postprob,postprob >= postprob.ordered[max(which(FDR.Bayes > FDR.threshold))])
      P3D_return_dat <- as.data.frame(P3D_return_dat)
      colnames(P3D_return_dat) <- c("Beta_Hat","Var_Beta_Hat","P_Values","PostProb","Significant")
    }

    if(sum(P3D_return_dat$Significant) > 0){
      indices_X <- which(P3D_return_dat$Significant)
      X <- SNPs[, indices_X,drop = FALSE]

      if(ite != 1){
        pi0=pi0_1
        g=g_1
      }else{
        pi0 = pi0.hat
        g = g.hat
      }


      PQL_est <- PQL(Y=Y, Z=Z, kinship=kinship, X=X, Xc=Xc, Xs=Xs, indices_X=indices_X, indices_Xc=indices_previous, family=family, postprob = postprob, offset =offset)
      y_star <- PQL_est$y_star
      Kappa <- PQL_est$kappa
      Beta <- PQL_est$beta
      inv_V_vector <- PQL_est$inv_v

      H <- PQL_est$H
      P <- PQL_est$P
      D_inv <- PQL_est$D_inv

      if(is.null(Xs)){
        y.tilde = t(P) %*% (y_star-Beta[1])

      }else{
        y.tilde = t(P) %*% (y_star-cbind(1,Xs)%*%Beta[1:(ncol(Xs)+1)])
      }
      X.tilde = t(P) %*% cbind(X, Xc)

      #y.tilda~N(X.tilde*Beta,D)  # full model, no intercept in X.tilde
      # elements with all k regressors
      dinvy <- D_inv*y.tilde
      ydinvy <- sum(D_inv*y.tilde^2)
      xdinvx <- t(X.tilde)%*%(D_inv*X.tilde)

      total.p <- ncol(X.tilde)

      if(total.p < 16){
        # Do full model search
        total.models <- 2^total.p
        log.unnormalized.posterior.probability <- rep(NA, total.models)
        log.unnormalized.posterior.probability[1] <- total.p * log(pi0) + log_marginal_likelihood_null(y.tilde=y.tilde,D_inv=D_inv)
        dat <- rep(list(0:1), total.p)
        dat <- as.matrix(expand.grid(dat))
        for (i in 1:(total.models-1)){
          model <- unname(which(dat[i + 1,] == 1))
          k <- length(model)

          Xsub <- X.tilde[,model,drop = FALSE]
          if(Matrix::rankMatrix(Xsub)[1] < ncol(Xsub)){
            dropped_cols <- caret::findLinearCombos(Xsub)$remove
            model <- model[-dropped_cols]
          }

          x.tilde_m <- matrix(X.tilde[,model], ncol = length(model))

          log.unnormalized.posterior.probability[i+1] <-
            k*log(1-pi0) + (total.p-k)*log(pi0) +
            log.marginal.likelihood(k=k, x.tilde_m=x.tilde_m, y.tilde=y.tilde, D_inv=D_inv, ydinvy=ydinvy, dinvy=dinvy, g=g)

        }
        log.unnormalized.posterior.probability <- log.unnormalized.posterior.probability - max(log.unnormalized.posterior.probability)
        unnormalized.posterior.probability <- exp(log.unnormalized.posterior.probability)
        posterior.probability <- unnormalized.posterior.probability/sum(unnormalized.posterior.probability)

      }else {
        # Do model search with genetic algorithm
        fitness_ftn <- function(string){
          if(sum(string) == 0){
            return(total.p * log(pi0) + log_marginal_likelihood_null(y.tilde=y.tilde,D_inv=D_inv))
          }else{
            model <- which(string==1)
            k <- length(model)

            Xsub <- X.tilde[,model,drop = FALSE]
            if(Matrix::rankMatrix(Xsub)[1] < ncol(Xsub)){
              dropped_cols <- caret::findLinearCombos(Xsub)$remove
              model <- model[-dropped_cols]
            }

            x.tilde_m <- matrix(X.tilde[,model], ncol = length(model))


            return(k*log(1-pi0) + (total.p-k)*log(pi0) +
                     log.marginal.likelihood(k=k, x.tilde_m=x.tilde_m, y.tilde=y.tilde, D_inv=D_inv, ydinvy=ydinvy, dinvy=dinvy, g=g)
            )
          }
        }

        if(total.p > 99){
          suggestedsol <- diag(total.p)
          tmp_log.unnormalized.posterior.probability <- vector()
          for(i in 1:total.p){
            model <- which(suggestedsol[i,]==1)
            k <- length(model)

            Xsub <- X.tilde[,model,drop = FALSE]
            if(Matrix::rankMatrix(Xsub)[1] < ncol(Xsub)){
              dropped_cols <- caret::findLinearCombos(Xsub)$remove
              model <- model[-dropped_cols]
            }

            x.tilde_m <- matrix(X.tilde[,model], ncol = length(model))

            tmp_log.unnormalized.posterior.probability[i] <- (k*log(1-pi0) + (total.p-k)*log(pi0) +
                                                                log.marginal.likelihood(k=k, x.tilde_m=x.tilde_m, y.tilde=y.tilde, D_inv=D_inv, ydinvy=ydinvy, dinvy=dinvy, g=g) )

          }
          suggestedsol <- rbind(0,suggestedsol[order(tmp_log.unnormalized.posterior.probability,decreasing = TRUE)[1:99],])
        }else{
          suggestedsol <- rbind(0,diag(total.p))
        }

        # maxiterations = 4000
        # runs_til_stop = 1000

        fitness_ftn <- memoise::memoise(fitness_ftn)
        ans <- GA::ga("binary", fitness = fitness_ftn, nBits = total.p,maxiter = maxiterations,popSize = 100,
                      elitism = min(c(10,2^total.p)),run = runs_til_stop,suggestions = suggestedsol,monitor = FALSE)
        memoise::forget(fitness_ftn)
        dat <- ans@population
        dupes <- duplicated(dat)
        dat <- dat[!dupes,]
        ans@fitness <- ans@fitness[!dupes]
        log.unnormalized.posterior.probability <- ans@fitness - max(ans@fitness)
        unnormalized.posterior.probability <- exp(log.unnormalized.posterior.probability)
        posterior.probability <- unnormalized.posterior.probability/sum(unnormalized.posterior.probability)
      }

      inclusion_prb <- unname((t(dat)%*%posterior.probability)/sum(posterior.probability))

      model <- dat[which.max(posterior.probability),]
      model_dat <- cbind(c(indices_X, indices_previous),model,inclusion_prb)
      model_dat <- as.data.frame(model_dat)
      colnames(model_dat) <- c("SNPs","BestModel","Inclusion_Prob")
      # select_list[[ite]] <- model_dat
      tmp=list(prescreen = P3D_return_dat,postprob=postprob,modelselection = model_dat,pi_0_hat = pi0.hat, g_hat = g.hat)

    }else{
      if(ite == 1){
        tmp=list(prescreen = P3D_return_dat,postprob=postprob,modelselection = "No significant in prescreen1",pi_0_hat = pi0.hat, g_hat = g.hat)
      }
    }


    if(!is.character(tmp$modelselection)){
      indices_Xc <- tmp$modelselection$SNPs[tmp$modelselection$BestModel == 1]
      postprob <- tmp$postprob
      indices_sig_list[[ite]] <- indices_Xc

      if(ite == 1){
        pi0_1 <- tmp$pi_0_hat
        g_1 <- tmp$g_hat
      }

      Xc <- SNPs[,indices_Xc]

      ite <- ite + 1
    }
  }

  return(tmp)

}

