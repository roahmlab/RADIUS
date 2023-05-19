#ifndef ROAHM_POINT_XY_HPP_
#define ROAHM_POINT_XY_HPP_

/// @file point_xy.hpp contains a struct to store (x, y) coordinates

namespace roahm {
/// Contains a \f$ (x, y) \f$ point, generally in meters
struct PointXY {
  /// the x position, generally in meters
  double x_;
  /// the y position, generally in meters
  double y_;

  /// Default constructor, defaults to \f$ (x, y) = (0\text{m}, 0\text{m}) \f$
  PointXY() : x_{}, y_{} {}

  /// Constructor given x and y position
  /// \param x the x position, generally in meters
  /// \param y the y position, generally in meters
  PointXY(double x, double y) : x_(x), y_(y) {}
};
} // namespace roahm
#endif // ROAHM_POINT_XY_HPP_
