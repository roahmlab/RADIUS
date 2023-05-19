#ifndef ROAHM_POINT_XYH_HPP_
#define ROAHM_POINT_XYH_HPP_
/// @file point_xyh.hpp contains a struct to store 2d poses (x, y, heading)
#include <ostream>

#include "simple_util.hpp"
namespace roahm {

/// Contains a \f$ (x, y, h) \f$ pose, generally in meters and radians
struct PointXYH {
  /// the x position, generally in meters
  double x_;
  /// the y position, generally in meters
  double y_;
  /// the heading, generally in radians
  double h_;
  /// Default constructor, defaults to
  /// \f$ (x, y, h) = (0\textrm{m}, 0\textrm{m}, 0\textrm{rad}) \f$
  constexpr PointXYH() : x_{}, y_{}, h_{} {}
  constexpr PointXYH(double x, double y, double h) : x_(x), y_(y), h_(h) {}

  /// Converts a point in the global frame into the local frame, mirroring
  /// about the y-axis if desired.
  /// \param new_frame the local frame pose w.r.t. the world frame
  /// \return the input pose converted from the world frame to the local frame
  [[nodiscard]] PointXYH ToLocalFrame(const PointXYH& new_frame) const {
    const double cos_h = std::cos(new_frame.h_);
    const double sin_h = std::sin(new_frame.h_);
    const double trans_x = x_ - new_frame.x_;
    const double trans_y = y_ - new_frame.y_;
    const double trans_rot_x = cos_h * trans_x - sin_h * trans_y;
    const double trans_rot_y = sin_h * trans_x + cos_h * trans_y;
    return {trans_rot_x, trans_rot_y, ToAngle(AngleDiff(h_, new_frame.h_))};
  }

  /// Mirrors the pose about the x-axis.
  /// \return the pose mirrored about the x-axis
  /// \f$ (x, y, h) \rightarrow (x, -y, -h) \f$
  [[nodiscard]] constexpr PointXYH Mirror() const { return {x_, -y_, -h_}; }

  friend std::ostream& operator<<(std::ostream& out, const PointXYH& pt) {
    return out << "[x: " << pt.x_ << "] [y: " << pt.y_ << "] [h: " << pt.h_
               << "]";
  }
};
} // namespace roahm
#endif // ROAHM_POINT_XYH_HPP_
