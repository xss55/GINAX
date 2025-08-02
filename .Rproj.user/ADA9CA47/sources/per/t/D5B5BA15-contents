
#' Performs GINA-X as described in the manuscript, Xu, Williams, Tegge, and Ferreira Genome-wide iterative fine-mapping for non-Gaussian data, Nature Genetics, Submitted.
#'
#'
#'
#' @param Y The observed phenotypes, count or binary.
#' @param SNPs The SNP matrix, where each column represents a single SNP encoded as the numeric coding 0, 1, 2. This is entered as a matrix object.
#' @param FDR_Nominal The nominal false discovery rate for which SNPs are selected from in the screening step.
#' @param family Specify if the response is count ("poisson") or binary ("bernoulli").
#' @param Covariance A list of covariance matrices that are the covariance matrices of the random effects. This matches the list of design matrices in Z.
#' @param Z A list of matrices specifying the design matrix of each random effect of interest.
#' @param offset If family = "poisson", the offset of each ecotype, can be a vector or a number if the number of offset is the same for each ecotype. If family = "binomial", offset = NULL.
#' @param maxiterations The maximum iterations the genetic algorithm in the model selection step iterates for, defaulted at 2000
#' @param runs_til_stop The number of iterations at the same best model before the genetic algorithm in the model selection step converges, defaulted at 400
#' @return The column indices of SNPs that were in the best model identified by GINAX
#' @examples
#' #library(GINAX)
#' #data("Y_binary");data("SNPs");data("kinship")
#' #n <- length(Y_binary)
#' #covariance <- list()
#' #covariance[[1]] <- kinship
#' #covariance[[2]] <- diag(1, nrow = n, ncol = n)
#'
#' #set.seed(1330)
#' #output_poisson <- GINAX(Y=Y_binary, SNPs=SNPs,
#'  #                  Covariance=covariance, Z=NULL, family="bernoulli",
#'   #                 offset=NULL, FDR_Nominal = 0.05,
#'    #                maxiterations = 2000, runs_til_stop = 400)
#'
#' @export
GINAX <- function(Y, Covariance, SNPs, family, Z=NULL, offset=NULL,
                  FDR_Nominal = 0.05, maxiterations = 2000, runs_til_stop = 400){

  #family <- match.arg(family)

  if(sum(family %in% c("poisson","bernoulli")) == 0){
    stop("family must be either poisson or bernoulli")
  }

  if(family == "bernoulli"){
    family <- "binomial"
  }

  if(!is.numeric(Y)){
    stop("Y has to be numeric")
  }
  if(!is.matrix(SNPs)){
    stop("SNPs has to be a matrix object")
  }
  if(!is.numeric(SNPs)){
    stop("SNPs has to contain numeric values")
  }
  if(maxiterations-floor(maxiterations)!=0){
    stop("maxiterations has to be a integer")
  }
  if(runs_til_stop-floor(runs_til_stop)!=0){
    stop("runs_til_stop has to be a integer")
  }
  if(maxiterations < runs_til_stop){
    stop("maxiterations has to be larger than runs_til_stop")
  }
  if(FDR_Nominal > 1 | FDR_Nominal < 0){
    stop("FDR_Nominal has to be between 0 and 1")
  }

  if(!is.null(offset)){
    offset <- exp(offset)
  }

  n <- length(Y)
  if(is.null(Z)){
    n_rf <- length(Covariance)
    Z <- list()
    for(rf in 1:n_rf){
      Z[[rf]] <- diag(1, ncol = n, nrow = n)
    }
  }

  GINAX <- GINAX_terminal(Y = Y, kinship = Covariance, Z=Z, SNPs=SNPs, family=family, offset=offset,
                             FDR.threshold = 1-FDR_Nominal, maxiterations = maxiterations, runs_til_stop = runs_til_stop)

  if(is.character(GINAX$modelselection)){
    return("No significant SNP")
  }else{
    return(sort(GINAX$modelselection$SNPs[GINAX$modelselection$BestModel == 1]))
  }

}
