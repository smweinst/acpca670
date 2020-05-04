#' @title Function to simulate data from the AC-PCA paper
#' @description Code for the simulation function based on AC-PCA paper
#' @name sim_dat_fun
#' @param n the number of subjects; default is 5
#' @param b the number of brain regions; default is 10
#' @param p the number of features per brain region; default is 400
#' @param alpha the constant that is multiplied by the confounder; default is 2.5
#' @return list
#' @export
sim_dat_fun = function(n=5,b=10,p=400,alpha=2.5){
  mu = 1:b
  w1 = (mu-mean(mu))/sd(mu)

  sigma_mat = matrix(nrow=b,ncol=b)

  for (i in 1:b){
    for (j in 1:b){
      sigma_mat[i,j] = exp(-((w1[i] - w1[j])^2)/4)
    }
  }

  w2 = MASS::mvrnorm(1,mu=rep(0,b),Sigma = 0.25*sigma_mat)
  W = cbind(w1,w2)

  # h is a 2 x p matrix whose rows are N(0,I_p)
  h = MASS::mvrnorm(2,mu=rep(0,p),Sigma = diag(1,p))
  Omega = W%*%h # dim: b x p

  # donor-specific componenent Gamma:
  ## G.i = L1.i + L2.i
  ### L1.i = 1r.i, 1 is vector of 1's, r.i is N(0,I_p)

  X.mat = matrix(nrow=0, ncol = p)
  Gamma.mat = matrix(nrow=0,ncol=p)
  for (i in 1:n){
    r.i = matrix(rnorm(p,0,1),nrow = 1, ncol = p)
    L1.i = matrix(rep(1,b),nrow=b) %*% r.i

    B.i = sample(c(rep(1,3), rep(0,b-3))) # B.i has 3 entries that are Uniform(0,2) and the rest are set to 0
    B.i[B.i==1] = runif(3,min=0,max=2) # sample from U(0,2) for the non-zero elements of B.i
    s.i = matrix(rnorm(p,0,1),nrow = 1, ncol = p)
    L2.i = B.i%*%s.i

    G.i = L1.i + L2.i

    X.i = Omega + alpha*G.i + rnorm(n=1,mean=0,sd=sqrt(0.25)) # epsilon (subject-level random noise)

    X.mat = rbind(X.mat,X.i)
    Gamma.mat = rbind(Gamma.mat,G.i)
  }

  labels = rep(1:b,n) # labels for each brain region (10 per person)

  group = rep(1:n,b)
  group = sort(group) # labels for each brain (i.e., each subject)

  # defing the confounder matrix, Y
  # this definition for Y is taken from code provided by author of AC-PCA paper
  Y <- c()
  for (k in 1:b){ # region
    for (i in 1:(n-1)){ # subject
      for (j in (i+1):n){ #
        tmp <- rep(0, n*b)
        tmp[(i-1)*b+k] <- 1
        tmp[(j-1)*b+k] <- -1
        Y <- cbind(Y, as.matrix(tmp))
      }
    }
  }

  # return a list with each part of the simulated data
  return(list(X.mat = X.mat,
              Gamma.mat=Gamma.mat,
              Y=Y,
              Omega = Omega,
              labels = labels,
              group=group))
}
