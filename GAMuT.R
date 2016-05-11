## GAMuT.R
##----------------------------------------------------------------------------------------------
## GAMuT example script
## Written 20160313
##
## Within this code, we define
## PHENO:    matrix of phenotypes (no missing data allowed)
## GENO:     matrix of rare variants (no missing data allowed)
## COVAR:    matrix of covariates (no missing data allowed)
##
## Simulated data included for 1000 unrelated subjects:
## 1) phenotypes.csv = matrix of 6 continuous phenotypes
## 2) variants.csv   = matrix of 353 rare variants
## 3) covariates.csv = matrix of 2 continuous covariates
##
##----------------------------------------------------------------------------------------------

##------------------------------------------------
## prerequisites for this code:
## * CompQuadForm
##   R package that provides davies() function,
##
## * GAMuT-functions.R
##   definitions for functions called in this code
##------------------------------------------------
library(CompQuadForm)
source("GAMuT-functions.R")


##===============================================================================
## STEP 1: read in PHENO, GENO, and COVAR matrices
##===============================================================================
PHENO <- read.csv("phenotypes.csv")
GENO  <- read.csv("variants.csv")
COVAR <- read.csv("covariates.csv")


##===============================================================================
## STEP 2:  residualize the phenotypes in PHENO on the covariates in COVAR
## Note: Each phenotype in PHENO will be residualized on each covariate in COVAR
##===============================================================================

##------------------------------------------------
## If you wish to only residualize phenotypes on a subset of covariates in COVAR, 
## make sure to make appropriate changes to the code below
##------------------------------------------------
residuals = matrix(NA, nrow=nrow(PHENO), ncol=ncol(PHENO)) # Matrix of residualized phenotypes
for(l in 1:ncol(PHENO)){
  model.residual = lm(PHENO[,l] ~ as.matrix(COVAR))
  residuals[,l] = resid(model.residual)
}


##===============================================================================
## STEP 3:  form the phenotypic similarity matrix Yc (using notation from the AJHG paper)
## and corresponding eigenvalues lambda_Y based on residualized phenotypes from Step 2
##===============================================================================
## centered & scaled matrix of residualized phenotypes:
P0 <- apply(residuals, 2, scale, scale=T, center=T)
P  <- as.matrix(P0) # convert dataframes into matrix

##------------------------------------------------
## for this example, we construct Yc using the projection matrix
##------------------------------------------------
## function for constructing the projection matrix and corresponding eigenvalues:
proj_pheno = proj_GAMuT_pheno(P)
Yc = proj_pheno$Kc                # Projection matrix
lambda_Y = proj_pheno$ev_Kc       # Eigenvalues of Yc

##------------------------------------------------
## If one wishes to construct Yc using the linear kernel instead of projection matrix, 
## use the following code instead 
##
## function for constructing the linear phenotype kernel and eigenvalues:
## linear_pheno <- linear_GAMuT_pheno(P) 
## Yc = linear_pheno$Kc           # Linear kernel similarity matrix
## lambda_Y = linear_pheno$ev_Kc  # Eigenvalues of Yc
##------------------------------------------------


##===============================================================================
## STEP 4:  form the genotypic similarity matrix Xc (using notation from the AJHG paper)
## and corresponding eigenvalues lambda_X
##===============================================================================

## We assume weighted linear kernel for genotypes where weights are function of 
## minor-allele frequency (MAF). We apply the beta-distribution MAF weights of 
## Wu et al. (2011) in the SKAT paper

## form the variant weights:
MAF = colMeans(GENO)/2                      # sample MAF of each variant in the sample
beta_weight = dbeta(MAF,1,25)/dbeta(0,1,25) # assume beta-distribution weights
## Note: one can use other weight functions by simply recoding beta_weight as desired

G0 = as.matrix(GENO)%*%diag(beta_weight) # Weighted rare variants
G  = as.matrix(scale(G0,center=T,scale=F)) # Centered genotype matrix

## function for constructing the weighted linear kernel matrix for genotypes and eigenvalues:
linear_geno <- linear_GAMuT_geno(G) 
Xc <- linear_geno$Lc            # Linear kernel similarity matrix
lambda_X <- linear_geno$ev_Lc   # Eigenvalues of Xc


##===============================================================================
## STEP 5: Construct the GAMuT Test and obtain the p-value
##===============================================================================
GAMuT_pvalue= TestGAMuT(Yc,lambda_Y,Xc,lambda_X)
