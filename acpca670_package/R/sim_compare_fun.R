#' @title Function to simulate data and compare different methods
#' @description Data will be simulated using sim_dat_fun function and then we will compare performance of ComBat, SVA, and acPCA, just like in the AC-PCA paper
#' @name sim_compare_fun
#' @param n the number of subjects; default is 5
#' @param b the number of brain regions; default is 10
#' @param p the number of features per brain region; default is 400
#' @param alpha the constant that is multiplied by the confounder; default is 2.5
#' @param nsim number of simulations to run; default is 100
#' @export
sim_compare_fun = function(n=5,b=10,p=400,alpha=2.5,nsim = 100){
  # initialize empty matrices to store correlations:
  pca.combat.scores.cor = pca.combat.loadings.cor = acpca.scores.cor = acpca.loadings.cor = pca.sva.scores.cor = pca.sva.loadings.cor = matrix(nrow=nsim,ncol=2)
  for (s in 1:nsim){
    sim_dat.s = sim_dat_fun(n=n,b=b,p=p,alpha=alpha)
    X.mat.s = sim_dat.s$X.mat
    omega.s = sim_dat.s$Omega
    omega.s.shared = do.call("rbind",replicate(n,omega.s,simplify = F))
    labels = sim_dat.s$labels
    group = sim_dat.s$group
    Y = sim_dat.s$Y
    # Gamma.mat.s = sim_dat.s$Gamma.mat

    # regular pca on the shared component:
    pca_omega = prcomp(omega.s.shared, center = T)
    pca_omega.scores = pca_omega$x # scores
    pca_omega.loadings = pca_omega$rotation # loadings

    # ComBat:
    mod.combat =model.matrix(~factor(labels))
    combat.X.s = sva::ComBat(t(X.mat.s),batch = group, mod = mod.combat) # transpose because sva package expects features in rows
    # apply pca after combat:
    pca.combat.s = prcomp(t(combat.X.s),center = T)
    pca.combat.scores.cor[s,] = sapply(1:2, FUN = function(t){
      cor(pca.combat.s$x[,t],pca_omega.scores[,t],method = "pearson") # correlation between scores
    })
    pca.combat.loadings.cor[s,] = sapply(1:2, FUN = function(t){
      cor(pca.combat.s$rotation[,t],pca_omega.loadings[,t], method = "spearman") # correlation between loadings
    })

    # SVA:
    sva.mod = model.matrix(~factor(labels))
    sva.X.s = sva::sva(t(X.mat.s), mod = sva.mod) # transpose because sva package expects features in rows
    sv = sva.X.s$sv
    fsva.X.s = sva::fsva(t(X.mat.s), sv = sva.X.s,mod = sva.mod,
                         newdat = t(X.mat.s))

    pca.sva.s = prcomp(t(fsva.X.s$db),center = T) # transpose again to get back to features in columns instead of rows

    pca.sva.scores.cor[s,] = sapply(1:2, FUN = function(t){
      cor(pca.sva.s$x[,t],pca_omega.scores[,t],method = "pearson") # correlation between scores
    })

    pca.sva.loadings.cor[s,] = sapply(1:2, FUN = function(t){
      cor(pca.sva.s$rotation[,t],pca_omega.loadings[,t],method = "pearson") # correlation between loadings
    })


    # AC-PCA:
    acpca.s.tune = acPCA::acPCAtuneLambda(X = X.mat.s,
                                          Y = Y,
                                          nPC = 2,
                                          lambdas = seq(0,10,0.05),
                                          anov=T, kernel = "linear",quiet = T)

    acpca.s = acPCA::acPCA(X = X.mat.s, Y = Y,
                           lambda = acpca.s.tune$best_lambda,
                           kernel = "linear", nPC = 2)

    acpca.scores.cor[s,] = sapply(1:2, FUN = function(t){
      cor(acpca.s$Xv[,t],pca_omega.scores[,t],method = "pearson") # correlation between scores
    })

    acpca.loadings.cor[s,] = sapply(1:2, FUN = function(t){
      cor(acpca.s$v[,t], pca_omega.loadings[,t],method = "spearman") # correlation between loadings
    })

  }

  if (nsim==1){
    par(mfrow=c(1,5))
    # true pattern:
    plot(pca_omega.scores[,1],pca_omega.scores[,2],type = "n", xlab = "PC 1", ylab = "PC 2",main = "True Pattern")
    text(pca_omega.scores[,1],pca_omega.scores[,2],labels = labels)

    # regular pca:
    reg_pca = prcomp(X.mat.s)
    plot(reg_pca$x[,1],reg_pca$x[,2],type = "n", xlab = "PC 1", ylab = "PC 2",main = "PCA")
    text(reg_pca$x[,1],reg_pca$x[,2],labels = labels,col = group+1)

    # combat:
    plot(pca.combat.s$x[,1],pca.combat.s$x[,2],type = "n",main = "ComBat",
         xlab = "PC1",ylab = "PC2")
    text(pca.combat.s$x[,1],pca.combat.s$x[,2],labels = labels,col = group+1)

    # sva:
    plot(pca.sva.s$x[,1],pca.sva.s$x[,2],type = "n",main = "SVA",xlab = "PC1",ylab = "PC2")
    text(pca.sva.s$x[,1],pca.sva.s$x[,2],labels = labels,col = group+1)

    # acpca:
    plot(acpca.s$Xv[,1],acpca.s$Xv[,2],type = "n",main = "ACPCA",xlab = "PC1",ylab = "PC2")
    text(acpca.s$Xv[,1],acpca.s$Xv[,2], labels = labels,col = group+1)

  } else {
    par(mfrow=c(1,1),mar=c(5,5,5,5),cex.axis=1,cex.main=1.5)
    # combat:
    vioplot::vioplot(
      # combat:
      abs(pca.combat.scores.cor[,1]),abs(pca.combat.loadings.cor[,1]),abs(pca.combat.scores.cor[,2]),abs(pca.combat.loadings.cor[,2]),

      # sva:
      abs(pca.sva.scores.cor[,1]),abs(pca.sva.loadings.cor[,1]), abs(pca.sva.scores.cor[,2]),abs(pca.sva.loadings.cor[,2]),

      # acpca:
      abs(acpca.scores.cor[,1]),abs(acpca.loadings.cor[,1]),abs(acpca.scores.cor[,2]),abs(acpca.loadings.cor[,2]),
      ylim = c(0,1),main = "",
      col = "white",
      names= rep(c("PC1","PC1","PC2","PC2"),3),
      cex.axis=1);mtext(
        side=1,at=1:12,line=1.75,text = rep(c("", "loading")))
    abline(v=c(4.5,8.5),lty=2)
    title("ComBat",line = 1,adj = 0.15)
    title("SVA",line = 1,adj = 0.5)
    title("AC-PCA",line = 1,adj = 0.85)
  }

}

