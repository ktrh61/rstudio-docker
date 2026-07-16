// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// X: N x K  （行=サンプル、列=共変量）※R側で x_w = sweep(design, 1, w_hat, `*` 済み
// W: G x N  （遺伝子×サンプルの重み、正）
// Y: G x N  （log2(tumor) - log2(normal)）
// sigma2: 長さ G （hatsg2）
// return_full_cov: 共分散のフル行列を返すか（TRUE ）/ 対角のみ（FALSE）

// [[Rcpp::export]]
Rcpp::List gls_per_gene_cpp(const arma::mat& X,
                            const arma::mat& W,
                            const arma::mat& Y,
                            const arma::vec& sigma2,
                            const bool return_full_cov = true) {
  const int N = X.n_rows;
  const int K = X.n_cols;
  const int G = Y.n_rows;
  
  if ((int)W.n_rows != G || (int)W.n_cols != N)
    Rcpp::stop("W shape mismatch");
  if ((int)Y.n_cols != N)
    Rcpp::stop("Y shape mismatch");
  if ((int)sigma2.n_elem != G)
    Rcpp::stop("sigma2 length mismatch");
  
  arma::mat coef(G, K, arma::fill::zeros);
  
  arma::cube cov_full;
  arma::mat  cov_diag;
  if (return_full_cov) {
    cov_full.set_size(K, K, G);
    cov_full.zeros();
  } else {
    cov_diag.set_size(G, K);
    cov_diag.zeros();
  }
  
  arma::mat XWX(K, K);
  arma::vec XWy(K);
  arma::mat Xw(N, K); // 行ごとに重み付けした X
  
  for (int g = 0; g < G; ++g) {
    arma::vec wg = arma::trans(W.row(g));   // N x 1
    arma::vec yg = arma::trans(Y.row(g));   // N x 1
    
    // 行ごとの重み付け：Xw = diag(wg) * X  をBLASに乗せる
    arma::mat Xw = X.each_col() % wg;       // N x K
    
    // XWX = X' * Xw
    arma::mat XWX = X.t() * Xw;             // K x K
    
    // XWy = X' * (wg .* y)
    arma::vec XWy = X.t() * (wg % yg);      // K x 1
    
    // Cholesky & solve
    arma::mat R;
    if (!arma::chol(R, XWX)) {
      Rcpp::stop(std::string("chol(XWX) failed at gene ") + std::to_string(g + 1));
    }
    arma::vec z = arma::solve(arma::trimatl(R.t()), XWy);
    arma::vec beta = arma::solve(arma::trimatu(R), z);
    coef.row(g) = beta.t();
    
    // 共分散：inv_sympd は potri 相当で速いことが多い
    if (return_full_cov) {
      arma::mat XWX_inv = arma::inv_sympd(XWX);
      cov_full.slice(g) = XWX_inv * sigma2[g];
    } else {
      arma::mat XWX_inv = arma::inv_sympd(XWX);
      for (int j = 0; j < K; ++j) cov_diag(g, j) = XWX_inv(j, j) * sigma2[g];
    }
  }
  
  if (return_full_cov) {
    return Rcpp::List::create(
      Rcpp::_["coef"] = coef,
      Rcpp::_["cov_full"] = cov_full
    );
  } else {
    return Rcpp::List::create(
      Rcpp::_["coef"] = coef,
      Rcpp::_["cov_diag"] = cov_diag
    );
  }
}
