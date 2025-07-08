#include "RcppArmadillo.h"
using namespace arma;

vec log_choose_vec(const ivec& n, const ivec& k) {
  if (n.n_elem != k.n_elem) {
    Rcpp::stop("Vectors 'n' and 'k' must have the same length.");
  }

  vec n_d = conv_to<vec>::from(n);
  vec k_d = conv_to<vec>::from(k);
  vec nk_d = n_d - k_d;

  vec result(n_d.n_elem, fill::none);

  // Logical masks
  uvec is_valid = (k >= 0) && (k <= n);
  uvec is_equal = (n == k);
  uvec is_zero  = (k == 0);
  uvec compute  = is_valid && is_equal == 0 && is_zero == 0;

  // Set results
  result.fill(-datum::inf);  // default for invalid cases

  // n == k → log(choose(n, n)) = 0
  result.elem(find(is_equal)).zeros();

  // k == 0 → log(choose(n, 0)) = 0
  result.elem(find(is_zero)).zeros();

  // General case
  uvec idx = find(compute);
  if (!idx.is_empty()) {
    vec n_sub  = n_d.elem(idx);
    vec k_sub  = k_d.elem(idx);
    vec nk_sub = nk_d.elem(idx);

    vec logc = lgamma(n_sub + 1.0) - lgamma(k_sub + 1.0) - lgamma(nk_sub + 1.0);
    result.elem(idx) = logc;
  }

  return result;
}

mat log_choose_mat(const imat& N, const imat& K) {
  if (size(N) != size(K)) {
    Rcpp::stop("N and K must be the same size.");
  }

  int n_rows = N.n_rows;
  int n_cols = N.n_cols;
  mat result(n_rows, n_cols, fill::none);

  for (int col = 0; col < n_cols; ++col) {
    ivec n_col = N.col(col);
    ivec k_col = K.col(col);

    // Use the optimized log_choose_vector function
    vec logc = log_choose_vec(n_col, k_col);

    result.col(col) = logc;
  }

  return result;
}

vec softmax(const vec &x) {
  vec y(x.size());
  y = exp(x);
  y(y.size() - 1) = 1; 
  y = y / (1.0 * sum(y));
  return y;
}

//' Compute softmax column-wise with fixed numerator = 1 for last entry of each segment
//' 
//' NOTE: Lambda should be returned so that the entry from second to last to last is always 1
//' 
//' @param G: M x N matrix (G = Phi %*% t(Psi))
//' @param rho: temperature parameter
//' @return: M x N matrix Lambda of softmax values
mat matrix_softmax(const mat& G, double rho) {
  // G has one row less than lambda - ie, G_S-1
  int M = G.n_rows;
  int N = G.n_cols;
  int S = (1.0 + sqrt(1.0 + (8.0 * (M+1.0)))) / 2.0;

  ivec segment_lengths(S-1);
  segment_lengths = reverse(regspace<ivec>(2, S-1));

  mat Lambda(M+1, N, fill::ones);
  int start = 0;

  for (int i = 0; i < S-2; ++i) {
    int len = segment_lengths(i);
    int end = start + len;

    for (int col = 0; col < N; ++col) {
      vec g_block = G.submat(start, col, end - 1, col);
      vec scaled = rho * g_block;

      vec exp_vals(len);
      exp_vals(span(0, len - 2)) = exp(scaled(span(0, len - 2)) - max(scaled));  // stabilization
      exp_vals(len - 1) = 1.0;

      double denom = 1.0 + accu(exp_vals);
      for (int k = 0; k < len; ++k) {
        Lambda(start + k, col) = exp_vals(k) / denom;
      }
    }

    start += len;
  }

  return Lambda;
}

uword od_size_to_stops(uword N) {
  return (1 + sqrt(1 + 8 * N)) / 2;
}
