# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

L1proj <- function(v, z, lambda0) {
    .Call(`_wavegrpreg_L1proj`, v, z, lambda0)
}

L12proj <- function(eta, M, p, z, lambda0) {
    .Call(`_wavegrpreg_L12proj`, eta, M, p, z, lambda0)
}

spgaug <- function(W, y, Sigma, sizes, z, theta, start, gamma, a_min, a_max, sig1, sig2, Msave, max1, max2, tol, L1proj) {
    .Call(`_wavegrpreg_spgaug`, W, y, Sigma, sizes, z, theta, start, gamma, a_min, a_max, sig1, sig2, Msave, max1, max2, tol, L1proj)
}

SPG_gbridge <- function(W, y, Sigma, sizes, z, lambda0, init, maxiter, tol0, gamma, a_min, a_max, sig1, sig2, Msave, max1, max2, tol, L1proj) {
    .Call(`_wavegrpreg_SPG_gbridge`, W, y, Sigma, sizes, z, lambda0, init, maxiter, tol0, gamma, a_min, a_max, sig1, sig2, Msave, max1, max2, tol, L1proj)
}

SPG_glasso <- function(W, y, Sigma, sizes, z, lambda0, init, gamma, a_min, a_max, sig1, sig2, Msave, max1, max2, tol, L12proj) {
    .Call(`_wavegrpreg_SPG_glasso`, W, y, Sigma, sizes, z, lambda0, init, gamma, a_min, a_max, sig1, sig2, Msave, max1, max2, tol, L12proj)
}

SPG_lasso <- function(W, y, Sigma, z, lambda0, init, gamma, a_min, a_max, sig1, sig2, Msave, max1, max2, tol, L1proj) {
    .Call(`_wavegrpreg_SPG_lasso`, W, y, Sigma, z, lambda0, init, gamma, a_min, a_max, sig1, sig2, Msave, max1, max2, tol, L1proj)
}

