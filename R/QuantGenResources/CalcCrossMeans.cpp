/*
 * Functions for calculating genomic estimates of
 * mean cross performance based on a GS model fitting
 * additive and dominance (digenic) effects. Effects
 * are modeled using AlphaSimR's coding scheme.
 * This version is OpenMP-parallelized to speed up diallel loops.
 */

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat calcCrossMean(arma::mat& geno,
                        arma::vec& a,
                        arma::vec& d,
                        arma::uword ploidy) {

  // Determine dimensions and allocate output
  arma::uword nInd     = geno.n_rows;
  arma::uword nSnp     = geno.n_cols;
  arma::uword nCrosses = nInd * (nInd - 1) / 2;
  arma::mat  output(nCrosses, 3, arma::fill::zeros);

  // Transpose genotype matrix and convert to integer type
  arma::umat genoT = arma::conv_to<arma::umat>::from(geno.t());

  // Precompute additive and dominance maps according to ploidy
  arma::vec x = arma::regspace<arma::vec>(0, ploidy) / double(ploidy);
  arma::vec xa = 2.0 * x - 1.0;
  arma::vec xd = -4.0 * (x % x) + 4.0 * x;

  arma::mat gam;
  if (ploidy == 2) {
    gam = { {2,0}, {1,1}, {0,2} };
    gam /= 2.0;
  } else if (ploidy == 4) {
    gam = { {6,0,0}, {3,3,0}, {1,4,1}, {0,3,3}, {0,0,6} };
    gam /= 6.0;
  } else if (ploidy == 6) {
    gam = {
      {20,0,0,0}, {10,10,0,0}, {4,12,4,0},
      {1,9,9,1}, {0,4,12,4}, {0,0,10,10}, {0,0,0,20}
    };
    gam /= 20.0;
  } else {
    Rcpp::stop("No gamete probabilities for this ploidy");
  }

  arma::mat mapA(ploidy+1, ploidy+1, arma::fill::zeros);
  arma::mat mapD(ploidy+1, ploidy+1, arma::fill::zeros);

  // Build lookup tables for parental genotype combinations
  for (arma::uword i = 0; i <= ploidy; ++i) {
    for (arma::uword j = 0; j <= ploidy; ++j) {
      arma::mat F = gam.row(i).t() * gam.row(j);
      for (arma::uword k = 0; k <= ploidy/2; ++k) {
        for (arma::uword l = 0; l <= ploidy/2; ++l) {
          mapA(i,j) += F(k,l) * xa(k+l);
          mapD(i,j) += F(k,l) * xd(k+l);
        }
      }
    }
  }

  // Parallel half-diallel cross computation
#pragma omp parallel for collapse(2) schedule(dynamic)
  for (arma::uword i = 0; i < nInd - 1; ++i) {
    for (arma::uword j = i + 1; j < nInd; ++j) {
      // Compute flat row index k for this pair
      arma::uword k = i * nInd - (i * (i + 1)) / 2 + (j - i - 1);
      double val = 0.0;
      for (arma::uword m = 0; m < nSnp; ++m) {
        val += mapA(genoT(m,i), genoT(m,j)) * a(m)
        + mapD(genoT(m,i), genoT(m,j)) * d(m);
      }
      output(k,0) = double(i + 1);
      output(k,1) = double(j + 1);
      output(k,2) = val;
    }
  }

  return output;
}
