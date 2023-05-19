#ifndef ROAHM_MU_SIGMA_MULTI_
#define ROAHM_MU_SIGMA_MULTI_

#include <cmath>
#include <limits>
#include <map>
#include <vector>

namespace roahm {
struct TimeMicroseconds {
public:
  long time_us_;

private:
  constexpr TimeMicroseconds() = delete;
  constexpr TimeMicroseconds(long time_us) : time_us_{time_us} {}

public:
  constexpr static TimeMicroseconds
  FromMicroseconds(const int64_t duration_us) {
    return TimeMicroseconds{duration_us};
  }

  template <typename T>
  constexpr static TimeMicroseconds FromMilliseconds(const T duration_ms)
    requires(std::is_same_v<T, int64_t>)
  {
    return TimeMicroseconds{duration_ms * 1000};
  }

  template <typename T>
  constexpr static TimeMicroseconds FromSeconds(const T duration_s)
    requires(std::is_same_v<T, double>)
  {
    return TimeMicroseconds{
        static_cast<long>(std::round(duration_s * 1'000'000))};
  }

  [[nodiscard]] constexpr double ToSecondsDouble() const {
    return static_cast<double>(time_us_) / 1'000'000.0;
  }

  constexpr static TimeMicroseconds Get10ms() { return {10000}; }

  [[nodiscard]] constexpr std::strong_ordering
  operator<=>(const TimeMicroseconds& oth) const {
    return time_us_ <=> oth.time_us_;
  }

  /*
  [[nodiscard]] constexpr bool operator<(const TimeMicroseconds& oth) const
      noexcept(false) {
    return time_us_ < oth.time_us_;
  }

  [[nodiscard]] constexpr bool operator==(const TimeMicroseconds& oth) const
      noexcept(false) {
    return time_us_ == oth.time_us_;
  }
  */

  TimeMicroseconds& operator+=(const TimeMicroseconds& oth) {
    time_us_ += oth.time_us_;
    return *this;
  }

  TimeMicroseconds operator*(const int& m) const { return {time_us_ * m}; }

  [[nodiscard]] int64_t operator/(const TimeMicroseconds& oth) const {
    return time_us_ / oth.time_us_;
  }
};

struct MuSigmaMulti {
private:
  const TimeMicroseconds us_between_mu_sigmas_;
  std::map<TimeMicroseconds, std::vector<double>> mu_sigma_mapping_;

public:
  constexpr static int kEltsPerMuSigma = 5;

  template <std::size_t N>
  MuSigmaMulti(const std::array<std::vector<double>, N>& mu_sigmas_incr,
               const TimeMicroseconds us_between_mu_sigmas)
      : us_between_mu_sigmas_{us_between_mu_sigmas}, mu_sigma_mapping_{} {
    for (int i = 0; i < N; ++i) {
      const TimeMicroseconds t{us_between_mu_sigmas_ * (i + 1)};
      mu_sigma_mapping_[t] = mu_sigmas_incr.at(i);
    }
  }

  // Returns the minimum number of mu sigma sets contained in any of the mu
  // sigma vectors Note that this is not the minimum number of numeric elements,
  // but rather sets of mu sigmas
  [[nodiscard]] std::size_t MinNumMuSigmas() const;

  const std::vector<double>& GetMuSigmaVec(const TimeMicroseconds& time) const;

  [[nodiscard]] TimeMicroseconds MicrosecondsBetweenMuSigmas() const;
};

} // namespace roahm
#endif // ROAHM_MU_SIGMA_MULTI_