// CDM_fast3_arma_enhanced.cpp
// 進捗表示付き・高速（RcppArmadillo版）
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List CDM_fast3_arma(const arma::mat& X, 
                          std::string random = "False", 
                          bool verbose = true) {
  
  const int d = X.n_rows;
  const int n = X.n_cols;
  
  if (verbose) {
    Rcpp::Rcout << "CDM_fast3_arma starting..." << std::endl;
    Rcpp::Rcout << "Data dimensions: " << d << " genes × " << n << " samples" << std::endl;
  }
  
  if (n <= 3 || d <= 1) Rcpp::stop("Invalid dimensions");
  
  const int n1 = (n + 1) / 2;
  const int n2 = n - n1;
  int r = std::min(n2 - 1, d);
  
  if (verbose) {
    Rcpp::Rcout << "Sample split: " << n1 << " + " << n2 << " samples" << std::endl;
  }
  
  // インデックス（原型の挙動に合わせる）
  arma::uvec index = arma::regspace<arma::uvec>(0, n - 1);
  if (random != "False") {
    if (verbose) Rcpp::Rcout << "Applying random split..." << std::endl;
    index = arma::shuffle(index);
  }
  arma::uvec idx1 = index.head(n1);
  arma::uvec idx2 = index.tail(n2);
  
  // 分割 & 行方向中心化（各特徴ごと）
  if (verbose) Rcpp::Rcout << "Splitting data and centering..." << std::endl;
  arma::mat X1 = X.cols(idx1);
  arma::mat X2 = X.cols(idx2);
  X1.each_col() -= arma::mean(X1, 1);
  X2.each_col() -= arma::mean(X2, 1);
  
  // 交差行列（正規化付き）
  if (verbose) Rcpp::Rcout << "Computing Cross-Data Matrix..." << std::endl;
  const double denom = std::sqrt((n1 - 1.0) * (n2 - 1.0));
  arma::mat Sd = (X1.t() * X2) / denom;
  
  // 必要ランクのみ
  r = std::min(r, (int)Sd.n_rows);
  r = std::min(r, (int)Sd.n_cols);
  
  if (verbose) {
    Rcpp::Rcout << "Cross-Data Matrix dimensions: " << Sd.n_rows << "×" << Sd.n_cols << std::endl;
    Rcpp::Rcout << "Computing " << r << " components via SVD..." << std::endl;
  }
  
  // SVD（必要成分のみ）
  arma::mat U, V;
  arma::vec s;
  arma::svd_econ(U, s, V, Sd, "both", "std");
  U = U.cols(0, r - 1);
  V = V.cols(0, r - 1);
  s = s.head(r);
  
  if (verbose) {
    Rcpp::Rcout << "SVD completed. Singular values range: " 
                << s.min() << " to " << s.max() << std::endl;
  }
  
  // 符号合わせ（原型の u^T X1^T X2 v >= 0 になるように）
  if (verbose) Rcpp::Rcout << "Applying sign correction..." << std::endl;
  arma::mat S = X1.t() * X2;
  for (int i = 0; i < r; ++i) {
    const double v = arma::dot(U.col(i), S * V.col(i));
    if (v < 0.0) V.col(i) *= -1.0;
  }
  
  // 固有ベクトル（特徴空間）再構成
  if (verbose) Rcpp::Rcout << "Reconstructing eigenvectors in feature space..." << std::endl;
  arma::vec inv_scale1 = 1.0 / arma::sqrt(s * (n1 - 1.0));
  arma::vec inv_scale2 = 1.0 / arma::sqrt(s * (n2 - 1.0));
  arma::mat H1 = X1 * (U * arma::diagmat(inv_scale1));
  arma::mat H2 = X2 * (V * arma::diagmat(inv_scale2));
  arma::mat VEC = 0.5 * (H1 + H2);
  
  // L2正規化
  arma::vec cnrm = arma::sqrt(arma::sum(arma::square(VEC), 0)).t();
  for (int j = 0; j < r; ++j) {
    if (cnrm(j) > 1e-10) VEC.col(j) /= cnrm(j);
  }
  
  // スコア（サンプル空間）
  if (verbose) Rcpp::Rcout << "Computing sample scores..." << std::endl;
  arma::mat Sc1 = U * arma::diagmat(arma::sqrt(n1 * s));
  arma::mat Sc2 = V * arma::diagmat(arma::sqrt(n2 * s));
  arma::mat scores_block(n, r, arma::fill::zeros);
  scores_block.rows(idx1) = Sc1;
  scores_block.rows(idx2) = Sc2;
  
  if (verbose) {
    Rcpp::Rcout << "CDM analysis completed successfully!" << std::endl;
    Rcpp::Rcout << "Components computed: " << r << std::endl;
    if (r > 0) {
      double var_explained = s(0) * s(0) / arma::sum(s % s) * 100;
      Rcpp::Rcout << "PC1 variance explained: " << std::fixed << std::setprecision(1) 
                  << var_explained << "%" << std::endl;
    }
  }
  
  // values は numeric ベクトルで返す（原型と型を厳密一致）
  return Rcpp::List::create(
    Rcpp::Named("values")  = Rcpp::NumericVector(s.begin(), s.end()),
    Rcpp::Named("vectors") = VEC,
    Rcpp::Named("scores")  = scores_block
  );
}