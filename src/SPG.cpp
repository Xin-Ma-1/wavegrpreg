#include <RcppArmadillo.h>
// [[Rcpp::depends( RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List SPG(arma::mat W, arma::vec y, arma::mat Sigma, arma::vec sizes, 
               double z, double lambda,
               double gamma, double a_min, double a_max, double sig1, 
               double sig2, int Msave, int max1, int max2, double tol){
  arma::vec start_idx = cumsum(sizes)-sizes;
  arma::vec end_idx   = cumsum(sizes)-1;
  int M               = start_idx.n_elem;
  int p               = W.n_cols;
  
  // initialization
  arma::vec eta = arma::zeros(M*p);
  arma::cube wtw = arma::zeros(p,p,M);
  arma::cube wty = arma::zeros(p,1,M);
  for(int m = 0; m < M; ++m){
    arma::mat w_m = W.rows(start_idx(m),end_idx(m));
    arma::vec y_m = y.subvec(start_idx(m),end_idx(m));
    wtw.slice(m) = w_m.t()*w_m/sizes(m) - Sigma;
    wty.slice(m) = w_m.t()*y_m/sizes(m);
  }
  // initial value of alpha
  arma::vec df = arma::zeros(M*p);
  for(int m = 0; m < M; ++m){
    arma::vec eta_m = eta.subvec(m*p, (m+1)*p-1);
    df.subvec(m*p, (m+1)*p-1) = wtw.slice(m)*eta_m - wty.slice(m);
  }
  arma::vec g1 = L12proj(eta-df, M, p, z, lambda) - eta;
  double alpha = 1/max(abs(g1));
  
  
  
  Rcpp::Rcout << "wtw:" << wtw << "\n" << "wty:" << wty  << "\n";
  
  
  
  
  return Rcpp::List::create(Rcpp::Named("start_idx")   = start_idx,
                            Rcpp::Named("end_idx")     = end_idx,
                            Rcpp::Named("M")           = M,
                            Rcpp::Named("p")           = p);
}