#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
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

// The two per-value predicates the counters above apply element-wise. The
// exact enumeration reuses them so that both paths agree on every boundary.
bool is_at_least_near(const double value, const double observed) {
  return value >= observed || near_equal(value, observed);
}

// Sorts one sample (group x first, `nx` of them) into ascending tie blocks and
// records which pooled positions belong to group x.
void tie_blocks_and_labels(const double* values, const int n, const int nx,
                           std::vector<int>& block_sizes,
                           std::vector<unsigned char>& in_x) {
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [values](const int lhs, const int rhs) {
              if (values[lhs] < values[rhs]) {
                return true;
              }
              if (values[rhs] < values[lhs]) {
                return false;
              }
              return lhs < rhs;
            });

  block_sizes.clear();
  in_x.assign(n, 0);
  int begin = 0;
  while (begin < n) {
    int end = begin + 1;
    while (end < n && values[order[end]] == values[order[begin]]) {
      ++end;
    }
    block_sizes.push_back(end - begin);
    for (int position = begin; position < end; ++position) {
      in_x[position] = static_cast<unsigned char>(order[position] < nx);
    }
    begin = end;
  }
}

// ---------------------------------------------------------------------------
// Exact enumeration over all C(n, nx) label allocations
// ---------------------------------------------------------------------------

// Binomial coefficient. Every value this is asked for during unranking is
// bounded by C(n, nx) itself, which callers have already capped, so unsigned
// 64-bit arithmetic cannot overflow here.
std::uint64_t binomial(const int n, int k) {
  if (k < 0 || k > n) {
    return 0;
  }
  k = std::min(k, n - k);
  std::uint64_t result = 1;
  for (int i = 1; i <= k; ++i) {
    result = result * static_cast<std::uint64_t>(n - k + i) /
             static_cast<std::uint64_t>(i);
  }
  return result;
}

// C(n, nx) as a double, or infinity once the count stops being exactly
// representable. Used to decide whether enumeration is affordable at all.
double n_allocations(const int n, const int nx) {
  const int k = std::min(nx, n - nx);
  double result = 1.0;
  for (int i = 1; i <= k; ++i) {
    result = result * static_cast<double>(n - k + i) / static_cast<double>(i);
    if (result > 9007199254740992.0) {  // 2^53
      return R_PosInf;
    }
  }
  return std::floor(result + 0.5);
}

// The `rank`-th combination of nx positions out of n, in the same lexicographic
// order that `advance_combination` below walks. Lets each worker thread start
// part-way through the enumeration.
void unrank_combination(const int n, const int nx, std::uint64_t rank,
                        std::vector<int>& out) {
  int lower = 0;
  for (int i = 0; i < nx; ++i) {
    const int last = n - (nx - i);
    int value = lower;
    for (; value < last; ++value) {
      const std::uint64_t count = binomial(n - 1 - value, nx - 1 - i);
      if (rank < count) {
        break;
      }
      rank -= count;
    }
    out[i] = value;
    lower = value + 1;
  }
}

bool advance_combination(const int n, const int nx, std::vector<int>& idx) {
  int i = nx - 1;
  while (i >= 0 && idx[i] == n - nx + i) {
    --i;
  }
  if (i < 0) {
    return false;
  }
  ++idx[i];
  for (int j = i + 1; j < nx; ++j) {
    idx[j] = idx[j - 1] + 1;
  }
  return true;
}

// Per-worker tallies. `bucket[c]` counts allocations for which exactly the
// first `c` sorted observed keys are extreme; a suffix sum turns the buckets
// into a per-gene exceedance count.
struct ExactTally {
  std::vector<std::int64_t> bucket;
  std::int64_t positive_infinite = 0;
  std::int64_t negative_infinite = 0;

  explicit ExactTally(const std::size_t n_keys) : bucket(n_keys + 1, 0) {}
};

// Walks allocation ranks [begin, end) and tallies them. Touches no R API, so
// it is safe to run on a worker thread.
void enumerate_range(const BMStatistic& statistic, const int n, const int nx,
                     const std::vector<double>& sorted_keys,
                     const int alternative, const double cran_scale,
                     const std::uint64_t begin, const std::uint64_t end,
                     ExactTally& tally) {
  const std::size_t n_keys = sorted_keys.size();
  std::vector<int> idx(nx);
  unrank_combination(n, nx, begin, idx);
  std::vector<unsigned char> in_x(n, 0);

  for (std::uint64_t r = begin; r < end; ++r) {
    std::fill(in_x.begin(), in_x.end(), 0);
    for (int i = 0; i < nx; ++i) {
      in_x[idx[i]] = 1;
    }
    const double raw = statistic(in_x) * cran_scale;
    if (raw == R_PosInf) {
      ++tally.positive_infinite;
    } else if (raw == R_NegInf) {
      ++tally.negative_infinite;
    }

    // Each alternative is folded into the single question "is `value` at least
    // as extreme as the key?": two-sided compares magnitudes, `greater` flips
    // the sign so that its `at most` test becomes an `at least` test.
    const double value = (alternative == 1)   ? std::fabs(raw)
                         : (alternative == 2) ? -raw
                                              : raw;

    std::size_t cut =
        static_cast<std::size_t>(std::lower_bound(sorted_keys.begin(),
                                                  sorted_keys.end(),
                                                  value + kCranTolerance) -
                                 sorted_keys.begin());
    // lower_bound works on a shifted threshold, which can land one position
    // off the element-wise predicate at the tolerance boundary. The predicate
    // is monotone in the key, so a local walk restores the exact cut.
    while (cut > 0 && !is_at_least_near(value, sorted_keys[cut - 1])) {
      --cut;
    }
    while (cut < n_keys && is_at_least_near(value, sorted_keys[cut])) {
      ++cut;
    }
    ++tally.bucket[cut];

    if (!advance_combination(n, nx, idx)) {
      break;
    }
  }
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

  std::vector<int> block_sizes;
  std::vector<unsigned char> observed_in_x;
  tie_blocks_and_labels(&data[0], n, nx, block_sizes, observed_in_x);

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

// [[Rcpp::export]]
double bm_n_allocations_cpp(const int n, const int nx) {
  if (nx < 0 || nx > n) {
    Rcpp::stop("nx must lie between zero and n");
  }
  return n_allocations(n, nx);
}

// Observed statistics and tie patterns for a gene-by-sample matrix, in one
// pass. The tie pattern is what determines the permutation null, so callers
// group genes by the returned pattern and evaluate each group's null once.
//
// [[Rcpp::export]]
Rcpp::List bm_observed_matrix_cpp(const Rcpp::NumericMatrix& data,
                                  const int nx) {
  const R_xlen_t n_long = data.ncol();
  if (n_long > std::numeric_limits<int>::max()) {
    Rcpp::stop("the combined sample is too large for the C++ implementation");
  }
  const int n = static_cast<int>(n_long);
  if (nx < 2 || nx > n - 2) {
    Rcpp::stop("each group must contain at least two observations");
  }
  const R_xlen_t n_genes = data.nrow();

  Rcpp::NumericVector statistics(n_genes);
  Rcpp::CharacterVector patterns(n_genes);
  std::vector<double> row(n);
  std::vector<int> block_sizes;
  std::vector<unsigned char> in_x;

  for (R_xlen_t gene = 0; gene < n_genes; ++gene) {
    for (int column = 0; column < n; ++column) {
      const double value = data(gene, column);
      if (!R_FINITE(value)) {
        Rcpp::stop("data must contain only finite values");
      }
      row[column] = value;
    }
    tie_blocks_and_labels(row.data(), n, nx, block_sizes, in_x);
    BMStatistic statistic(block_sizes, nx);
    statistics[gene] = statistic(in_x);

    std::string pattern;
    pattern.reserve(block_sizes.size() * 3);
    for (std::size_t i = 0; i < block_sizes.size(); ++i) {
      if (i > 0) {
        pattern.push_back(',');
      }
      pattern += std::to_string(block_sizes[i]);
    }
    patterns[gene] = pattern;
  }

  return Rcpp::List::create(Rcpp::Named("statistic") = statistics,
                            Rcpp::Named("pattern") = patterns);
}

// Exact permutation counts for every supplied observed statistic at once: a
// single enumeration of all C(n, nx) allocations answers all of them, because
// the permutation null depends on the data only through the tie pattern.
//
// [[Rcpp::export]]
Rcpp::List bm_exact_extreme_cpp(const Rcpp::IntegerVector& block_sizes_r,
                                const int nx,
                                const Rcpp::NumericVector& observed_standard,
                                const int alternative,
                                const int n_threads) {
  if (alternative < 1 || alternative > 3) {
    Rcpp::stop("invalid alternative code");
  }
  std::vector<int> block_sizes = Rcpp::as<std::vector<int>>(block_sizes_r);
  BMStatistic statistic(block_sizes, nx);
  const int n = std::accumulate(block_sizes.begin(), block_sizes.end(), 0);
  const int ny = n - nx;
  if (nx < 2 || ny < 2) {
    Rcpp::stop("each group must contain at least two observations");
  }

  const double total_double = n_allocations(n, nx);
  if (!R_FINITE(total_double)) {
    Rcpp::stop("the number of allocations is too large to enumerate exactly");
  }
  const std::uint64_t total = static_cast<std::uint64_t>(total_double);

  const double cran_scale =
      static_cast<double>(n) /
      (static_cast<double>(nx) * static_cast<double>(ny));

  // Fold the alternative into a single "key" per gene, so that a gene is
  // extreme exactly when some allocation's folded value is at least its key.
  const R_xlen_t n_genes = observed_standard.size();
  std::vector<double> keys(n_genes);
  for (R_xlen_t gene = 0; gene < n_genes; ++gene) {
    const double observed = observed_standard[gene] * cran_scale;
    if (ISNAN(observed)) {
      Rcpp::stop("observed statistic must not be NaN");
    }
    keys[gene] = (alternative == 1)   ? std::fabs(observed)
                 : (alternative == 2) ? -observed
                                      : observed;
  }

  // Infinite keys cannot be resolved by threshold arithmetic, so they are held
  // back and answered from the infinity tallies once enumeration finishes.
  std::vector<R_xlen_t> finite_genes;
  finite_genes.reserve(n_genes);
  for (R_xlen_t gene = 0; gene < n_genes; ++gene) {
    if (R_FINITE(keys[gene])) {
      finite_genes.push_back(gene);
    }
  }
  std::sort(finite_genes.begin(), finite_genes.end(),
            [&keys](const R_xlen_t lhs, const R_xlen_t rhs) {
              return keys[lhs] < keys[rhs];
            });
  std::vector<double> sorted_keys(finite_genes.size());
  for (std::size_t i = 0; i < finite_genes.size(); ++i) {
    sorted_keys[i] = keys[finite_genes[i]];
  }

  int workers = n_threads;
  if (workers < 1) {
    workers = 1;
  }
  // Threading only pays off once each worker has real work to do.
  const std::uint64_t minimum_chunk = 200000;
  if (total / static_cast<std::uint64_t>(workers) < minimum_chunk) {
    workers = static_cast<int>(
        std::max<std::uint64_t>(1, total / minimum_chunk));
  }

  std::vector<ExactTally> tallies(workers, ExactTally(sorted_keys.size()));
  if (workers == 1) {
    enumerate_range(statistic, n, nx, sorted_keys, alternative, cran_scale, 0,
                    total, tallies[0]);
  } else {
    std::vector<std::thread> threads;
    threads.reserve(workers);
    const std::uint64_t chunk = total / static_cast<std::uint64_t>(workers);
    for (int worker = 0; worker < workers; ++worker) {
      const std::uint64_t begin = chunk * static_cast<std::uint64_t>(worker);
      const std::uint64_t end =
          (worker == workers - 1) ? total
                                  : chunk * static_cast<std::uint64_t>(worker + 1);
      threads.emplace_back([&, begin, end, worker]() {
        enumerate_range(statistic, n, nx, sorted_keys, alternative, cran_scale,
                        begin, end, tallies[worker]);
      });
    }
    for (auto& thread : threads) {
      thread.join();
    }
  }

  std::int64_t positive_infinite = 0;
  std::int64_t negative_infinite = 0;
  std::vector<std::int64_t> bucket(sorted_keys.size() + 1, 0);
  for (const ExactTally& tally : tallies) {
    positive_infinite += tally.positive_infinite;
    negative_infinite += tally.negative_infinite;
    for (std::size_t i = 0; i < bucket.size(); ++i) {
      bucket[i] += tally.bucket[i];
    }
  }

  // A gene at sorted position j is extreme for every allocation whose cut
  // exceeded j, so the per-gene counts are a suffix sum over the buckets.
  Rcpp::NumericVector counts(n_genes);
  std::int64_t running = 0;
  for (std::size_t position = sorted_keys.size(); position >= 1; --position) {
    running += bucket[position];
    counts[finite_genes[position - 1]] = static_cast<double>(running);
  }

  const double value_positive_infinite =
      (alternative == 1)
          ? static_cast<double>(positive_infinite + negative_infinite)
          : (alternative == 2 ? static_cast<double>(negative_infinite)
                              : static_cast<double>(positive_infinite));
  for (R_xlen_t gene = 0; gene < n_genes; ++gene) {
    if (R_FINITE(keys[gene])) {
      continue;
    }
    // A key of -Inf is met by every allocation; a key of +Inf only by those
    // allocations whose own folded value is infinite.
    counts[gene] = (keys[gene] == R_PosInf) ? value_positive_infinite
                                            : total_double;
  }

  return Rcpp::List::create(Rcpp::Named("count") = counts,
                            Rcpp::Named("total") = total_double);
}
