#ifndef ROAHM_DYN_OBS_HPP_
#define ROAHM_DYN_OBS_HPP_
#include <ostream>

#include "point_xyh.hpp"

namespace roahm {
class DynObs {
public:
  enum class ObstacleType {
    kStaticBoundary, // We won't use any mu sigma stuff here
    kDynamicObs
  };

private:
  double c_x0_;
  double c_y0_;
  double heading_rad_;
  double velocity_;
  double length_;
  double width_;
  ObstacleType obs_type_;

public:
  double CenterX0() const { return c_x0_; }
  double CenterY0() const { return c_y0_; }
  double HeadingRad() const { return heading_rad_; }
  double Velocity() const { return velocity_; }
  double Length() const { return length_; }
  double Width() const { return width_; }
  DynObs(double x, double y, double h, double v, double l, double w,
         ObstacleType obs_type);
  DynObs RelativeTo(PointXYH pt, bool mirror) const;
  [[nodiscard]] bool IsStaticBoundary() const;

  friend std::ostream& operator<<(std::ostream& out, const DynObs& obs);
  bool operator==(const DynObs& other) const;
};

} // namespace roahm
#endif // ROAHM_DYN_OBS_HPP_