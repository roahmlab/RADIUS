#include "mu_sigma_multi.hpp"
namespace roahm {
[[nodiscard]] std::size_t MuSigmaMulti::MinNumMuSigmas() const {
  if (mu_sigma_mapping_.empty()) {
    return 0;
  }
  std::size_t ret_val = std::numeric_limits<decltype(ret_val)>::max();
  for (const auto& [key, val] : mu_sigma_mapping_) {
    ret_val = std::min(val.size(), ret_val);
  }
  return ret_val / kEltsPerMuSigma;
}

const std::vector<double>&
MuSigmaMulti::GetMuSigmaVec(const TimeMicroseconds& time) const {
  return mu_sigma_mapping_.at(time);
}

[[nodiscard]] TimeMicroseconds
MuSigmaMulti::MicrosecondsBetweenMuSigmas() const {
  return us_between_mu_sigmas_;
}
} // namespace roahm