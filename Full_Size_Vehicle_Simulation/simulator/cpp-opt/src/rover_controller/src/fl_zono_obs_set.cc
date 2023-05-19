#include "fl_zono_obs_set.hpp"
#include <cassert>
#include <fmt/format.h>
namespace roahm {
FlZonoObsSet::FlZonoObsSet() : dyn_obs_{} {}
FlZonoObsSet FlZonoObsSet::RelativeTo(PointXYH pt, bool mirror) const {
  FlZonoObsSet ret;
  for (const auto& obs : dyn_obs_) {
    ret.PushObs(obs.RelativeTo(pt, mirror));
  }
  return ret;
}

void FlZonoObsSet::ClearAll() { dyn_obs_.clear(); }
std::size_t FlZonoObsSet::GetNumObs() const { return dyn_obs_.size(); }
void FlZonoObsSet::PushObs(DynObs obs) { dyn_obs_.push_back(obs); }
void FlZonoObsSet::PrintObsConv(const int obs_idx,
                                const ::roahm::Interval t_interval) const {
  const auto obs = dyn_obs_.at(obs_idx);
  const auto center = ComputeCenter(obs_idx, t_interval);
  const auto gens = GetGenerators(obs_idx, t_interval);
  std::cout << "Obs [" << obs_idx << "] at t = [" << t_interval.Min() << ", "
            << t_interval.Max() << "]" << std::endl;
  std::cout << "| Obs Params: " << obs << std::endl;
  std::cout << "| Center: " << center.x_ << ", " << center.y_ << std::endl;
  std::cout << "| Gen 1: " << gens.x1_ << ", " << gens.y1_ << std::endl;
  std::cout << "| Gen 2: " << gens.x2_ << ", " << gens.y2_ << std::endl;
}
FlZonoObsSet::CenterXY
FlZonoObsSet::ComputeCenter(const int obs_idx,
                            const ::roahm::Interval t_interval) const {
  // TODO move this into GetCenter
  const auto dyn_obs = dyn_obs_.at(obs_idx);
  const double t_center = t_interval.Midpoint();
  const double velocity = dyn_obs.Velocity();
  const double heading_rad = dyn_obs.HeadingRad();

  const double dist_from_start = velocity * t_center;
  const double cos_h = std::cos(heading_rad);
  const double sin_h = std::sin(heading_rad);
  const double x0 = dyn_obs.CenterX0();
  const double y0 = dyn_obs.CenterY0();
  return {x0 + dist_from_start * cos_h, y0 + dist_from_start * sin_h};
}

FlZonoObsSet::GenXY
FlZonoObsSet::GetGenerators(int obs_idx, ::roahm::Interval t_interval) const {
  const auto dyn_obs = dyn_obs_.at(obs_idx);
  const double length = dyn_obs.Length();
  const double velocity = dyn_obs.Velocity();
  const double width = dyn_obs.Width();
  const double heading_rad = dyn_obs.HeadingRad();

  const double delta_t = t_interval.Width();

  const double length_with_vel_adjusted{length + (velocity * delta_t)};
  const double longitudinal_gen_mag =
      (dyn_obs.IsStaticBoundary())
          ? (length_with_vel_adjusted / 2.0)
          : (((5.0 / 3.0) * length_with_vel_adjusted / 2.0 + 2.4));

  // The maximum extent (max - min) that an objects footprint could achieve
  // along the width dimension sampling at 3 sigma.
  constexpr double kFootprintWidthExtent{2.2 + 10.0 * (1.32 / 3)};
  // Footprint gen (width dim): width / 2.0
  // Sampling gen (width dim): kFootprintWidthExtent / 2.0 - width / 2.0
  // Total gen: Footprint + Sampling == kFootprintWidthExtent / 2.0-
  const double lateral_gen_mag = (dyn_obs.IsStaticBoundary())
                                     ? (width / 2.0)
                                     : (kFootprintWidthExtent / 2.0);

  // fmt::print("Width: {} (FPW: {}) (Y0: {}) (long G: {}) (lat G: {}) (length:
  // {}) (vel: {}) (delta_t: {}) (IsStatic: {})\n",
  //   width, kFootprintWidthExtent, dyn_obs.CenterY0(), longitudinal_gen_mag,
  //   lateral_gen_mag, length, velocity, delta_t, dyn_obs.IsStaticBoundary());
  // const double lateral_gen_mag = (width / 2.0);
  const double cos_h = std::cos(heading_rad);
  const double sin_h = std::sin(heading_rad);
  const double x_gen_1 = longitudinal_gen_mag * cos_h;
  const double y_gen_1 = longitudinal_gen_mag * sin_h;
  const double x_gen_2 = lateral_gen_mag * -sin_h;
  const double y_gen_2 = lateral_gen_mag * cos_h;
  return {x_gen_1, y_gen_1, x_gen_2, y_gen_2};
}

FlZonoObsSet::CenterXY
FlZonoObsSet::GetSingleGenerator(int obs_idx, int gen_idx,
                                 ::roahm::Interval t_interval) const {
  assert(gen_idx >= 0);
  assert(gen_idx < 2);
  assert(obs_idx < static_cast<int>(GetNumObs()));
  auto gens = GetGenerators(obs_idx, t_interval);
  if (gen_idx == 0) {
    return {gens.x1_, gens.y1_};
  }
  return {gens.x2_, gens.y2_};
}
[[nodiscard]] DynObs FlZonoObsSet::GetSingleDynObs(int obs_idx) const {
  return dyn_obs_.at(obs_idx);
}
} // namespace roahm
