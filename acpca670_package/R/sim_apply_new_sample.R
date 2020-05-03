#' @title Function to apply loadings from AC-PCA to new sample
#' @description User can specify alpha in the original sample and alpha in the new sample (i.e., samples differ in terms of strength of confounder)
#' @name sim_apply_new_sample
#' @param n the number of subjects; default is 5
#' @param b the number of brain regions; default is 10
#' @param p the number of features per brain region; default is 400
#' @param alpha_orig the constant that is multiplied by the confounder in the original simulated sample; default is 2.5
#' @param alpha_new the constant that is multiplied by the confounder in the new simulated sample; default is also 2.5
#' @param nsim number of simulations to run; default is 100
#' @export
sim_apply_new_sample = function(n=5,b=10,p=400,alpha_orig=2.5,alpha_new=2.5,nsim = 100){
  scores.cor.omega.new = matrix(nrow=nsim,ncol=2)
  for (s in 1:nsim){
    # simulate original dataset from which AC-PCA loadings will be obtained
    sim_dat.s = sim_dat_fun(n=n,b=b,p=p,alpha=alpha_orig)
    X.mat.s = sim_dat.s$X.mat
    Y = sim_dat.s$Y

    # tune lambda
    acpca.s.tune = acPCA::acPCAtuneLambda(X = X.mat.s,
                                          Y = Y,
                                          nPC = 2,
                                          lambdas = seq(0,10,0.05),
                                          anov=T, kernel = "linear",quiet = T)
    # get acpca loadings
    acpca.s.loadings = acPCA::acPCA(X = X.mat.s, Y = Y,
                                    lambda = acpca.s.tune$best_lambda,
                                    kernel = "linear", nPC = 2)$v[,1:2]

    # simulate a new dataset with alpha_new
    sim_dat.s.new = sim_dat_fun(n=n,b=b,p=p,alpha=alpha_new)
    X.mat.s.new = sim_dat.s.new$X.mat

    # apply AC-PCA loadings from original sample to new sample:
    Xv.newsamp = X.mat.s.new%*%acpca.s.loadings

    omega.s.new = sim_dat.s.new$Omega
    omega.s.new.shared = do.call("rbind",replicate(n,omega.s.new,simplify = F))

    # pca on omega in new data:
    pca_omega.new.scores = prcomp(omega.s.new.shared, center = T)$x

    scores.cor.omega.new[s,] = sapply(1:2, FUN = function(t){
      cor(Xv.newsamp[,t], pca_omega.new.scores[,t],method = "pearson") # correlation between scores
    })

  }
  if (nsim > 1){
    par(mfrow=c(1,1))
    vioplot::vioplot(abs(scores.cor.omega.new[,1]),abs(scores.cor.omega.new[,2]),
                     ylim = c(0,1), ylab = c("Pearson correlation"),
                     col = "white", names = c("PC1","PC2"),
                     main = "Correlation with shared component when AC-PCA from another sample is used");mtext(
                       bquote(paste(alpha['original']," = ",.(alpha_orig), "   ", alpha['new'], " = ", .(alpha_new))),side = 3
                       )
  }
  else{
    par(mfrow=c(1,2))
    plot(pca_omega.new.scores[,1],pca_omega.new.scores[,2], main = "True Pattern",
         xlab = "PC1", ylab = "PC2",type = 'n');text(
           pca_omega.new.scores[,1],pca_omega.new.scores[,2],labels = sim_dat.s.new$labels
         )
    plot(Xv.newsamp[,1],Xv.newsamp[,2], type = 'n',xlab = "PC1",ylab = "PC2",
         main = "AC-PCA from different sample applied");text(
           labels = sim_dat.s.new$labels, col = sim_dat.s.new$group + 1,
           Xv.newsamp[,1],Xv.newsamp[,2],
         );    mtext(
           bquote(paste(alpha['original']," = ",.(alpha_orig), "   ", alpha['new'], " = ", .(alpha_new))),side = 3
         )

  }

}
