#include "dyn_obs.hpp"

#include "gencon.hpp"

namespace roahm {
DynObs::DynObs(double x, double y, double h, double v, double l, double w,
               ::roahm::DynObs::ObstacleType obs_type)
    : c_x0_{x}, c_y0_{y}, heading_rad_{::roahm::ToAngle(h)}, velocity_{v},
      length_{l}, width_{w}, obs_type_{obs_type} {}

DynObs DynObs::RelativeTo(PointXYH pt, bool mirror) const {
  const double mirror_mult = mirror ? -1.0 : 1.0;
  const double dx = c_x0_ - pt.x_;
  const double dy = c_y0_ - pt.y_;
  const double dh = ::roahm::ToAngle(heading_rad_ - pt.h_);
  const double cos_h = std::cos(pt.h_);
  const double sin_h = std::sin(pt.h_);
  const double new_x = cos_h * dx + sin_h * dy;
  const double new_y = -sin_h * dx + cos_h * dy;
  DynObs ret{new_x,
             mirror_mult * new_y,
             mirror_mult * dh,
             velocity_,
             length_,
             width_,
             obs_type_};
  // std::cout << "Tranforming DynObs relative to " << pt << "[mirror: " <<
  // mirror << "]\n   " << *this << "\n-> " << ret << std::endl;
  return ret;
}

[[nodiscard]] bool DynObs::IsStaticBoundary() const {
  return obs_type_ == ObstacleType::kStaticBoundary;
}

std::ostream& operator<<(std::ostream& out, const DynObs& obs) {
  return std::cout << "[c_x0: " << obs.c_x0_ << "] [c_y0: " << obs.c_y0_
                   << "] [h: " << obs.heading_rad_
                   << "] [vel: " << obs.velocity_
                   << "] [length: " << obs.length_ << "] [width: " << obs.width_
                   << "]";
}

bool DynObs::operator==(const DynObs& other) const {
  return (c_x0_ == other.c_x0_) && (c_y0_ == other.c_y0_) &&
         (heading_rad_ == other.heading_rad_) &&
         (velocity_ == other.velocity_) && (length_ == other.length_) &&
         (width_ == other.width_);
}
} // namespace roahm