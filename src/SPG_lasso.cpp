#include <RcppArmadillo.h>
// [[Rcpp::depends( RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List SPG_lasso(arma::mat W, arma::vec y, arma::mat Sigma,
                     double z, double lambda0, arma::vec init,
                     double gamma, double a_min, double a_max, double sig1,
                     double sig2, int Msave, int max1, int max2, double tol,
                     Rcpp::Function L1proj){

  // initialization
  int p = W.n_cols;
  int n = y.n_elem;
  arma::vec eta = arma::zeros(p);
  if(is_finite(init)) eta = init;
  arma::mat wtw = W.t()*W/n - Sigma;
  arma::mat wty = W.t()*y/n;

  // initial value of alpha
  arma::vec df = wtw*eta - wty;
  Rcpp::NumericVector etadf = Rcpp::wrap(eta-df);
  Rcpp::NumericVector proj = L1proj(etadf, z, lambda0);
  arma::vec arma_proj = Rcpp::as<arma::vec>(proj);
  arma::vec g1 = arma_proj - eta;
  double g1max = max(abs(g1));
  if(g1max==0)Rcpp::Rcout << "Warning: constraint too strong!" << "\n";
  double alpha = 1/g1max;

  // main part of algorithm
  int iter = 0;
  int evaltime = 0;
  arma::vec fsave = arma::ones(Msave) * (-1e30);
  int fcount = 0;
  while(iter < max1 && evaltime < max2){
    // check if stationary
    arma::vec df = wtw*eta - wty;
    Rcpp::NumericVector etadf = Rcpp::wrap(eta-df);
    Rcpp::NumericVector proj = L1proj(etadf, z, lambda0);
    arma::vec arma_proj = Rcpp::as<arma::vec>(proj);
    g1 = arma_proj - eta;
    g1max = max(abs(g1));
    if(g1max<=tol)break;

    // backtracking
    double lambda= 1;
    Rcpp::NumericVector etadf1 = Rcpp::wrap(eta-alpha*df);
    Rcpp::NumericVector proj1 = L1proj(etadf1, z, lambda0);
    arma::vec arma_proj1 = Rcpp::as<arma::vec>(proj1);
    arma::vec dk = arma_proj1 - eta;
    arma::mat feta = 0.5*(eta.t()*wtw*eta) - wty.t()*eta;
    if(fcount < Msave){
      fsave(fcount) = feta(0,0);
      fcount++;
    }
    else{
      arma::vec fold = fsave;
      for(int i = 0; i < Msave-1; ++i){
        fsave(i) = fold(i+1);
      }
      fsave(Msave-1) = feta(0,0);
    }
    double fmax = max(fsave);
    double dkdf = sum(dk % df);
    arma::vec sk = arma::zeros(p);
    arma::vec yk = arma::zeros(p);
    // while loop to search for new eta
    while(1){
      arma::vec eta_plus = eta + lambda*dk;
      arma::mat feval = 0.5*(eta_plus.t()*wtw*eta_plus) - wty.t()*eta_plus;
      evaltime++;
      if(feval(0,0)<=fmax + gamma*lambda*dkdf){
        arma::vec eta_old = eta;
        eta = eta_plus;
        sk = eta - eta_old;
        arma::vec df_old = df;
        df = wtw*eta - wty;
        yk = df - df_old;
        break;
      }
      else{
        arma::vec randnum = arma::randu(1);
        lambda = (randnum(0,0)*(sig2-sig1) + sig1)*lambda;
      }
    }

    // update alpha
    double bk = sum(sk % yk);
    if(bk<=0) alpha = a_max;
    else{
      double ak = sum(sk % sk);
      alpha = ak/bk;
      if(alpha<a_min) alpha = a_min;
      if(alpha>a_max) alpha = a_max;
    }

    iter++;
  }

  return Rcpp::List::create(Rcpp::Named("eta")         = eta,
                            Rcpp::Named("z")           = z,
                            Rcpp::Named("lambda0")     = lambda0,
                            Rcpp::Named("iteration")   = iter,
                            Rcpp::Named("evaltime")    = evaltime,
                            Rcpp::Named("diff")        = g1max);
}
