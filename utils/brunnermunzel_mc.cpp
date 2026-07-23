#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace {

struct RunningVariance {
  int n = 0;
  double mean = 0.0;
  double m2 = 0.0;

  void add(const double value) {
    ++n;
    const double delta = value - mean;
    mean += delta / static_cast<double>(n);
    const double delta2 = value - mean;
    m2 += delta * delta2;
  }

  double sample_variance() const {
    if (n < 2) {
      return 0.0;
    }
    // Round-off can make a mathematically zero M2 very slightly negative.
    return std::max(0.0, m2 / static_cast<double>(n - 1));
  }
};

class BMStatistic {
 public:
  BMStatistic(const std::vector<int>& block_sizes, const int nx)
      : n_(std::accumulate(block_sizes.begin(), block_sizes.end(), 0)),
        nx_(nx),
        ny_(n_ - nx) {
    if (nx_ < 2 || ny_ < 2) {
      Rcpp::stop("each group must contain at least two observations");
    }
    int begin = 0;
    for (const int block_size : block_sizes) {
      if (block_size <= 0) {
        Rcpp::stop("tie block sizes must be positive integers");
      }
      const int end = begin + block_size;
      tie_blocks_.emplace_back(begin, end);
      begin = end;
    }
  }

  double operator()(const std::vector<unsigned char>& in_x) const {
    double pooled_sum_x = 0.0;
    double pooled_sum_y = 0.0;
    RunningVariance d_x;
    RunningVariance d_y;
    int seen_x = 0;
    int seen_y = 0;

    // The values do not change under label permutation. Therefore the sorted
    // tie blocks are cached, while pooled and within-group average ranks and
    // both variance components are calculated anew for every allocation.
    for (const auto& block : tie_blocks_) {
      const int begin = block.first;
      const int end = block.second;
      int block_x = 0;
      for (int pos = begin; pos < end; ++pos) {
        block_x += static_cast<int>(in_x[pos]);
      }
      const int block_y = end - begin - block_x;

      // begin is zero-based, end is one past the block. The first and last
      // one-based pooled ranks are begin + 1 and end.
      const double pooled_rank =
          (static_cast<double>(begin + 1) + static_cast<double>(end)) / 2.0;
      const double within_rank_x =
          static_cast<double>(seen_x) +
          (static_cast<double>(block_x) + 1.0) / 2.0;
      const double within_rank_y =
          static_cast<double>(seen_y) +
          (static_cast<double>(block_y) + 1.0) / 2.0;

      for (int pos = begin; pos < end; ++pos) {
        if (in_x[pos]) {
          pooled_sum_x += pooled_rank;
          d_x.add(pooled_rank - within_rank_x);
        } else {
          pooled_sum_y += pooled_rank;
          d_y.add(pooled_rank - within_rank_y);
        }
      }
      seen_x += block_x;
      seen_y += block_y;
    }

    const double mean_rank_x = pooled_sum_x / static_cast<double>(nx_);
    const double mean_rank_y = pooled_sum_y / static_cast<double>(ny_);
    const double s_x = d_x.sample_variance();
    const double s_y = d_y.sample_variance();
    const double numerator =
        static_cast<double>(nx_) * static_cast<double>(ny_) *
        (mean_rank_y - mean_rank_x);
    const double denominator =
        static_cast<double>(n_) *
        std::sqrt(static_cast<double>(nx_) * s_x +
                  static_cast<double>(ny_) * s_y);

    if (denominator > 0.0) {
      return numerator / denominator;
    }
    if (numerator > 0.0) {
      return R_PosInf;
    }
    if (numerator < 0.0) {
      return R_NegInf;
    }
    // In particular, all-identical observations map 0 / 0 to the neutral
    // statistic so every identical permutation is counted as extreme.
    return 0.0;
  }

 private:
  int n_;
  int nx_;
  int ny_;
  std::vector<std::pair<int, int>> tie_blocks_;
};

constexpr double kCranTolerance = 1e-14;

bool near_equal(const double lhs, const double rhs) {
  if (lhs == rhs) {
    return true;
  }
  if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
    return false;
  }
  return std::fabs(lhs - rhs) < kCranTolerance;
}

double count_at_least_near(const Rcpp::NumericVector& sorted,
                           const double observed) {
  const auto first_extreme = std::partition_point(
      sorted.begin(), sorted.end(), [observed](const double value) {
        return !(value >= observed || near_equal(value, observed));
      });
  return static_cast<double>(std::distance(first_extreme, sorted.end()));
}

double count_at_most_near(const Rcpp::NumericVector& sorted,
                          const double observed) {
  const auto first_non_extreme = std::partition_point(
      sorted.begin(), sorted.end(), [observed](const double value) {
        return value <= observed || near_equal(value, observed);
      });
  return static_cast<double>(
      std::distance(sorted.begin(), first_non_extreme));
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List bm_observed_cpp(const Rcpp::NumericVector& data, const int nx) {
  const R_xlen_t n_long = data.size();
  if (n_long > std::numeric_limits<int>::max()) {
    Rcpp::stop("the combined sample is too large for the C++ implementation");
  }
  const int n = static_cast<int>(n_long);
  if (nx < 2 || nx > n - 2) {
    Rcpp::stop("each group must contain at least two observations");
  }
  for (int i = 0; i < n; ++i) {
    if (!R_FINITE(data[i])) {
      Rcpp::stop("data must contain only finite values");
    }
  }

  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&data](const int lhs, const int rhs) {
    if (data[lhs] < data[rhs]) {
      return true;
    }
    if (data[rhs] < data[lhs]) {
      return false;
    }
    return lhs < rhs;
  });

  std::vector<int> block_sizes;
  std::vector<unsigned char> observed_in_x(n, 0);
  int begin = 0;
  while (begin < n) {
    int end = begin + 1;
    while (end < n && data[order[end]] == data[order[begin]]) {
      ++end;
    }
    block_sizes.push_back(end - begin);
    for (int position = begin; position < end; ++position) {
      observed_in_x[position] =
          static_cast<unsigned char>(order[position] < nx);
    }
    begin = end;
  }

  BMStatistic statistic(block_sizes, nx);
  return Rcpp::List::create(
      Rcpp::Named("statistic") = statistic(observed_in_x),
      Rcpp::Named("block_sizes") = Rcpp::wrap(block_sizes));
}

// [[Rcpp::export]]
Rcpp::NumericVector bm_mc_null_cpp(const Rcpp::IntegerVector& block_sizes_r,
                                   const int nx,
                                   const int B) {
  Rcpp::RNGScope rng_scope;

  if (B <= 0) {
    Rcpp::stop("B must be a positive integer");
  }
  std::vector<int> block_sizes = Rcpp::as<std::vector<int>>(block_sizes_r);
  BMStatistic statistic(block_sizes, nx);
  const int n = std::accumulate(block_sizes.begin(), block_sizes.end(), 0);
  const int ny = n - nx;
  const double cran_scale =
      static_cast<double>(n) /
      (static_cast<double>(nx) * static_cast<double>(ny));

  Rcpp::NumericVector null_distribution(B);
  std::vector<unsigned char> in_x(n, 0);
  for (int repetition = 0; repetition < B; ++repetition) {
    int remaining_x = nx;
    for (int position = 0; position < n; ++position) {
      const int remaining = n - position;
      bool choose_x = false;
      if (remaining_x == remaining) {
        choose_x = true;
      } else if (remaining_x > 0) {
        const double probability =
            static_cast<double>(remaining_x) / static_cast<double>(remaining);
        choose_x = R::runif(0.0, 1.0) < probability;
      }
      in_x[position] = static_cast<unsigned char>(choose_x);
      if (choose_x) {
        --remaining_x;
      }
    }

    null_distribution[repetition] = statistic(in_x) * cran_scale;
  }
  std::sort(null_distribution.begin(), null_distribution.end());
  return null_distribution;
}

// [[Rcpp::export]]
double bm_mc_extreme_cpp(const Rcpp::NumericVector& sorted_null,
                         const double observed_standard,
                         const int nx,
                         const int ny,
                         const int alternative) {
  if (nx < 2 || ny < 2) {
    Rcpp::stop("each group must contain at least two observations");
  }
  if (alternative < 1 || alternative > 3) {
    Rcpp::stop("invalid alternative code");
  }
  if (std::isnan(observed_standard)) {
    Rcpp::stop("observed statistic must not be NaN");
  }
  if (sorted_null.size() == 0) {
    Rcpp::stop("null distribution must not be empty");
  }

  // The returned observed statistic remains on the standard BM scale. CRAN's
  // permutation implementation uses this positive constant rescaling only for
  // its extreme comparison, which leaves all exact orderings unchanged.
  const double cran_scale =
      static_cast<double>(nx + ny) /
      (static_cast<double>(nx) * static_cast<double>(ny));
  const double observed = observed_standard * cran_scale;

  if (alternative == 2) {
    return count_at_most_near(sorted_null, observed);
  }
  if (alternative == 3) {
    return count_at_least_near(sorted_null, observed);
  }

  const double absolute_observed = std::fabs(observed);
  if (absolute_observed < kCranTolerance) {
    return static_cast<double>(sorted_null.size());
  }
  return count_at_most_near(sorted_null, -absolute_observed) +
         count_at_least_near(sorted_null, absolute_observed);
}
