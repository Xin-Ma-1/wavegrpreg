#include <RcppArmadillo.h>
// [[Rcpp::depends( RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec spgaug(arma::mat W, arma::vec y, arma::mat Sigma, arma::vec sizes,
                 double z, arma::vec theta, arma::vec start,
                 double gamma, double a_min, double a_max, double sig1,
                 double sig2, int Msave, int max1, int max2, double tol,
                 Rcpp::Function L1proj){
  arma::vec start_idx  = cumsum(sizes)-sizes;
  arma::vec end_idx    = cumsum(sizes)-1;
  int M                = start_idx.n_elem;
  int p                = W.n_cols;

  // initialization
  arma::vec eta = arma::zeros(M*p);
  if(is_finite(start)) eta = start;
  arma::cube wtw = arma::zeros(p,p,M);
  arma::cube wty = arma::zeros(p,1,M);
  for(int m = 0; m < M; ++m){
    arma::mat w_m = W.rows(start_idx(m),end_idx(m));
    arma::vec y_m = y.subvec(start_idx(m),end_idx(m));
    wtw.slice(m) = w_m.t()*w_m/sizes(m) - Sigma;
    for(int i=0; i<p; ++i) wtw.slice(m).col(i) = theta % wtw.slice(m).col(i);
    for(int j=0; j<p; ++j) wtw.slice(m).row(j) = theta.t() % wtw.slice(m).row(j);
    wty.slice(m) = (theta % (w_m.t()*y_m))/sizes(m);
  }

  // initial value of alpha
  arma::vec df = arma::zeros(M*p);
  for(int m = 0; m < M; ++m){
    arma::vec eta_m = eta.subvec(m*p, (m+1)*p-1);
    df.subvec(m*p, (m+1)*p-1) = wtw.slice(m)*eta_m - wty.slice(m);
  }
  Rcpp::NumericVector etadf = Rcpp::wrap(eta-df);
  Rcpp::NumericVector proj = L1proj(etadf, 1e30, 1);
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
    arma::vec df = arma::zeros(M*p);
    for(int m = 0; m < M; ++m){
      arma::vec eta_m = eta.subvec(m*p, (m+1)*p-1);
      df.subvec(m*p, (m+1)*p-1) = wtw.slice(m)*eta_m - wty.slice(m);
    }
    Rcpp::NumericVector etadf = Rcpp::wrap(eta-df);
    Rcpp::NumericVector proj = L1proj(etadf, 1e30, 1);
    arma::vec arma_proj = Rcpp::as<arma::vec>(proj);
    g1 = arma_proj - eta;
    g1max = max(abs(g1));
    if(g1max<=tol)break;

    // backtracking
    double lambda= 1;
    Rcpp::NumericVector etadf1 = Rcpp::wrap(eta-alpha*df);
    Rcpp::NumericVector proj1 = L1proj(etadf1, 1e30, 1);
    arma::vec arma_proj1 = Rcpp::as<arma::vec>(proj1);
    arma::vec dk = arma_proj1 - eta;
    arma::mat feta = arma::zeros(1,1);
    for(int m = 0; m < M; ++m){
      arma::vec eta_m = eta.subvec(m*p, (m+1)*p-1);
      feta = feta + 0.5*(eta_m.t()*wtw.slice(m)*eta_m) - wty.slice(m).t()*eta_m;
    }
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
    arma::vec sk = arma::zeros(M*p);
    arma::vec yk = arma::zeros(M*p);
    // while loop to search for new eta
    while(1){
      arma::vec eta_plus = eta + lambda*dk;
      arma::mat feval = arma::zeros(1,1);
      for(int m = 0; m < M; ++m){
        arma::vec eta_m = eta_plus.subvec(m*p, (m+1)*p-1);
        feval = feval + 0.5*(eta_m.t()*wtw.slice(m)*eta_m) - wty.slice(m).t()*eta_m;
      }
      evaltime++;
      if(feval(0,0)<=fmax + gamma*lambda*dkdf){
        arma::vec eta_old = eta;
        eta = eta_plus;
        sk = eta - eta_old;
        arma::vec df_old = df;
        df = arma::zeros(M*p);
        for(int m = 0; m < M; ++m){
          arma::vec eta_m = eta.subvec(m*p, (m+1)*p-1);
          df.subvec(m*p, (m+1)*p-1) = wtw.slice(m)*eta_m - wty.slice(m);
        }
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

  //from above we have eta.star, now get back to eta with L1 constraint at z
  arma::vec eta0 = arma::zeros(M*p);
  for(int m = 0; m < M; ++m){
    eta0.subvec(m*p, (m+1)*p-1) = theta % eta.subvec(m*p, (m+1)*p-1);
  }
  Rcpp::NumericVector etainput = Rcpp::wrap(eta0);
  Rcpp::NumericVector etaproj = L1proj(etainput, z, 0);
  eta = Rcpp::as<arma::vec>(etaproj);

  return eta;
}


// function to perform optimization with group bridge penalty
// [[Rcpp::depends( RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List SPG_gbridge(arma::mat W, arma::vec y, arma::mat Sigma, arma::vec sizes,
                       double z, double lambda0, arma::vec init,
                       int maxiter, double tol0,
                       double gamma, double a_min, double a_max, double sig1,
                       double sig2, int Msave, int max1, int max2, double tol,
                       Rcpp::Function L1proj){
  int M = sizes.n_elem;
  int p = W.n_cols;

  // initialization
  arma::vec eta = arma::zeros(M*p);
  if(is_finite(init)) eta = init;

  int iter = 0;
  double diff = 1;
  while(iter<=maxiter){
    arma::mat mat_eta(eta);
    mat_eta.reshape(p,M);

    // update augmented parameter theta
    arma::vec theta = lambda0*sqrt(sum(abs(mat_eta),1));

    // update eta with spgaug
    arma::vec eta_new = spgaug(W,y,Sigma,sizes,z,theta,eta,
                               gamma,a_min,a_max,sig1,sig2,Msave,max1,max2,tol,
                               L1proj);
    diff = max(abs(eta_new - eta));
    eta = eta_new;
    if(diff<tol0) break;
    iter++;
  }

  // convergence evaluation
  int iteration = 0;
  int convergence = 0;
  if(iter==1){
    iteration = iter;
  }else{
    if(iter>maxiter) iteration = iter-1;
    else{
      iteration = iter;
      convergence = 1;
    }
  }

  return Rcpp::List::create(Rcpp::Named("eta")         = eta,
                            Rcpp::Named("z")           = z,
                            Rcpp::Named("lambda0")     = lambda0,
                            Rcpp::Named("iteration")   = iteration,
                            Rcpp::Named("convergence") = convergence,
                            Rcpp::Named("diff")        = diff);
}
