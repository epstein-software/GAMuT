## GAMuT-functions.R
##
## 2016-March-13
##
## this script contains functions used in the main program
## for GAMuT analysis
##
##----------------------------------------------------------------------
## descriptions of individual functions:
##----------------------------------------------------------------------
##
## * proj_GAMuT_pheno
##   constructs projection matrix and corresponding eigenvalues for phenotypes
##
## * linear_GAMuT_pheno
##   constructs linear kernel matrix and corresponding eigenvalues for phenotypes
##
## * linear_GAMuT_geno
##   constructs linear kernel matrix and corresponding eigenvalues for genotypes 
##
## * TestGAMuT
##   constructs GAMuT statistic and returns p-value


##----------------------------------------------------------------------
## phenotypic similarity:  projection matrix
##----------------------------------------------------------------------
proj_GAMuT_pheno <- function(X){
    out = vector("list", 2)
    names(out) = c("Kc", "ev_Kc")

    out$Kc = X %*% solve(t(X) %*% X) %*% t(X)   # projection matrix
    out$ev_Kc = rep(1, ncol(X))                 # find eigenvalues
    return(out)	
}

##----------------------------------------------------------------------
## phenotypic similarity:  linear kernel
##----------------------------------------------------------------------
linear_GAMuT_pheno <- function(X){
    out = vector("list", 2)
    names(out) = c("Kc", "ev_Kc")
    
    Kc = t(X) %*% X   # transposed kernel to find eigenvalues
    ev_Kc = eigen(Kc, symmetric=T, only.values=T)$values  
    
    out$Kc = X %*% t(X) # similiarity kernel for test
    out$ev_Kc = ev_Kc[ev_Kc > 1e-08]
    return(out)
}

##----------------------------------------------------------------------
## genotypic similarity:  linear kernel
##----------------------------------------------------------------------
linear_GAMuT_geno <- function(X){
    out = vector("list", 2)
    names(out) = c("Lc", "ev_Lc")
    
    Lc = t(X) %*% X   # transposed kernel to find eigenvalues
    ev_Lc = eigen(Lc, symmetric=T, only.values=T)$values  
    
    out$Lc = X %*% t(X) # similiarity kernel for test
    out$ev_Lc = ev_Lc[ev_Lc > 1e-08]
    return(out)	 
}


##----------------------------------------------------------------------
## constructing the GAMuT statistic and deriving p-value:
##----------------------------------------------------------------------
TestGAMuT <- function(Yc, lambda_Y, Xc, lambda_X) {

    ## test statistic:
    m = nrow(Yc) # number of subjects in study
    GAMuT = (1/m) * sum(sum(t(Yc) * Xc))  
    

    ## populate vector of all pairwise combination of eigenvalues
    ## from the phenotype and genotype similarity matrices:
    Z <- (as.matrix(lambda_Y)) %*% t(as.matrix(lambda_X))
    Zsort <- sort(Z, decreasing=T)

    ## derive p-value of GAMuT statistic:
    scoredavies = GAMuT*m^2
    results_score <- davies(scoredavies, Zsort)
    davies_pvalue <- (results_score$Qq)

    return(davies_pvalue)
} 
