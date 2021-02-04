#include <Rcpp.h>
using namespace Rcpp;

// sign function
double sign(double x){
  if(x>0){
    return 1;
  } else if(x==0){
    return 0;
  } else{
    return -1;
  }
}

// cpp version of L1 projection with constraint
// [[Rcpp::export]]
NumericVector L1proj(NumericVector v, double z, double lambda0){
  int p = v.size();
  NumericVector abs_v = abs(v);
  NumericVector sign_v (p);
  for(int i = 0; i < p; ++i) sign_v[i] = sign(v[i]);
  NumericVector v_proj (p);
  NumericVector vec_zero (p);

  if(sum(abs_v)<=z){
    v_proj = pmax((abs_v-lambda0),vec_zero);
  }
  else{
    IntegerVector U = seq(0,(p-1));
    double s = 0;
    int rho = 0;
    while(U.size() > 0){
      int k = sample(U,1)[0];
      NumericVector vu = abs_v[U];
      double vk = abs_v[k];
      IntegerVector G = U[vu >= vk];
      IntegerVector L = setdiff(U, G);
      int drho = G.size();
      NumericVector vG = abs_v[G];
      double ds = sum(vG);
      if(s+ds-(rho+drho)*vk < z){
        s += ds;
        rho += drho;
        U = L;
      }
      else{
        IntegerVector setk = {k};
        U = setdiff(G, setk);
      }
    }
    double theta = (s-z)/rho;
    double lambda = lambda0;
    if(theta>lambda0) lambda = theta;
    v_proj = pmax((abs_v-lambda), vec_zero);
  }
  v_proj = v_proj*sign_v;

  return v_proj;
}


// cpp version of L-1,2 projection with constraint
// [[Rcpp::export]]
NumericVector L12proj(NumericVector eta, int M, int p, double z, double lambda0) {
  NumericMatrix mat_eta(p, M, eta.begin());
  NumericVector cnorm (p);
  for(int i=0; i<p; ++i){
    cnorm[i] = sqrt(sum(pow(mat_eta(i,_),2)));
  }
  NumericVector w = L1proj(cnorm, z, lambda0);
  NumericMatrix mat_eta_proj(p, M);
  for(int i=0; i<p; ++i){
    if(cnorm[i]>0)mat_eta_proj(i,_) = mat_eta(i,_)*w[i]/cnorm[i];
  }
  NumericVector eta_proj = as<NumericVector>(mat_eta_proj);
  eta_proj.attr("dim") = R_NilValue;

  return eta_proj;
}
