// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

/*
 * Calculate cross means (and segregation variance) for different ploidy levels
 * Output columns:
 *  1: parent i (1-based), 2: parent j (1-based), 3: mean, 4: within-cross variance
 */
// [[Rcpp::export]]
arma::mat calcCrossMean(arma::mat& geno,
                        arma::vec& a,
                        arma::vec& d,
                        arma::uword ploidy){
  arma::uword nInd = geno.n_rows;
  arma::uword nSnp = geno.n_cols;
  arma::uword nCrosses = nInd * (nInd-1) / 2;
  arma::mat output(nCrosses, 4, arma::fill::zeros);

  arma::umat genoT = arma::conv_to<arma::umat>::from(geno.t());

  // AlphaSimR genotype codings
  arma::vec x  = arma::regspace(0, ploidy);
  x /= double(ploidy);
  arma::vec xa = 2.0 * x - 1.0;                 // additive
  arma::vec xd = (-4.0) * (x % x) + 4.0 * x;    // digenic dominance  <-- THIS LINE

  // Gamete probabilities
  arma::mat gam;
  if(ploidy==2){
    gam = {{2,0},{1,1},{0,2}};
    gam /= 2.0;
  }else if(ploidy==4){
    gam = {{6,0,0},{3,3,0},{1,4,1},{0,3,3},{0,0,6}};
    gam /= 6.0;
  }else if(ploidy==6){
    gam = {{20,0,0,0},{10,10,0,0},{4,12,4,0},{1,9,9,1},{0,4,12,4},{0,0,10,10},{0,0,0,20}};
    gam /= 20.0;
  }else{
    Rcpp::stop("No gamete probabilities for this ploidy");
  }

  // Maps for E[xa], E[xd], and second moments per parental genotype pair
  arma::mat mapA(ploidy+1, ploidy+1, arma::fill::zeros);
  arma::mat mapD(ploidy+1, ploidy+1, arma::fill::zeros);
  arma::mat mapA2(ploidy+1, ploidy+1, arma::fill::zeros);
  arma::mat mapD2(ploidy+1, ploidy+1, arma::fill::zeros);
  arma::mat mapAD(ploidy+1, ploidy+1, arma::fill::zeros);

  for(arma::uword i=0; i<=ploidy; i++){
    for(arma::uword j=0; j<=ploidy; j++){
      arma::mat F = gam.row(i).t() * gam.row(j);
      for(arma::uword k=0; k<=ploidy/2; k++){
        for(arma::uword l=0; l<=ploidy/2; l++){
          arma::uword g = k + l;
          double w    = F(k,l);
          double xa_g = xa(g);
          double xd_g = xd(g);
          mapA(i,j)  += w * xa_g;
          mapD(i,j)  += w * xd_g;
          mapA2(i,j) += w * xa_g * xa_g;
          mapD2(i,j) += w * xd_g * xd_g;
          mapAD(i,j) += w * xa_g * xd_g;
        }
      }
    }
  }

  // Accumulate mean and within-cross variance across loci
  arma::uword k=0;
  for(arma::uword i=0; i<(nInd-1); i++){
    for(arma::uword j=i+1; j<nInd; j++){
      output(k,0) = i;
      output(k,1) = j;

      double mean_sum = 0.0;
      double var_sum  = 0.0;

      for(arma::uword m=0; m<nSnp; m++){
        arma::uword gi = genoT(m,i);
        arma::uword gj = genoT(m,j);

        double EA  = mapA(gi, gj);
        double ED  = mapD(gi, gj);
        double EA2 = mapA2(gi, gj);
        double ED2 = mapD2(gi, gj);
        double EAD = mapAD(gi, gj);

        double am = a(m);
        double dm = d(m);

        double mu_m  = EA*am + ED*dm;
        double ET2_m = (am*am)*EA2 + (dm*dm)*ED2 + 2.0*am*dm*EAD;
        double var_m = ET2_m - mu_m*mu_m;

        mean_sum += mu_m;
        var_sum  += var_m;
      }

      output(k,2) = mean_sum;
      output(k,3) = var_sum;
      k++;
    }
  }

  // 1-based parent indices for R
  output.col(0) += 1.0;
  output.col(1) += 1.0;

  return output;
}
