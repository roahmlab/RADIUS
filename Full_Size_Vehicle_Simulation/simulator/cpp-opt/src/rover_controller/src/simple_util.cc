#include "simple_util.hpp"

namespace roahm {

double Interval::DistanceTo(double val) const {
  // If val < min: val - min_ < 0
  // If val > max: max_ - val < 0
  const double max_dist = val - max_;
  const double min_dist = min_ - val;
  return std::max(std::max(max_dist, min_dist), 0.0);
}

} // namespace roahm