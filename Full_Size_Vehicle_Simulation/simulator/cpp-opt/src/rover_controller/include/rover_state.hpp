#ifndef ROAHM_COMPLEX_HEADING_STATE_HPP_
#define ROAHM_COMPLEX_HEADING_STATE_HPP_

#include <string> // for string

#include "point_xyh.hpp"
#include "simple_util.hpp"

/// @file rover_state.hpp Contains a class to store the rover's state
/// and some additional helper functions

namespace roahm {

/// Stores heading information as complex number in its
/// sine and cosine components to allow for averaging
/// multiple states, and guarantees that output heading
/// will remain within \f$ [-\pi, \pi] \f$
class RoverState {
public:
  /// x-coordinate in the global frame, in meters
  double x_;
  /// y-coordinate in the global frame, in meters
  double y_;

private:
  /// Heading, \f$ \theta \in [-\pi, \pi] \f$
  double h_;

public:
  /// Longitudinal velocity, forward/backwards in the robot frame [m/s]
  double u_;

  /// Lateral velocity, left/right in the robot frame [m]
  double v_;

  /// Rotational velocity, \f$\frac{\mathrm{d}}{\mathrm{d}t} \theta\f$ in the
  /// robot frame, [rads/s]
  double r_;

  /// Wheel speed, [m/s]
  double w_;

  /// Time-integrated angular velocity error from controller
  /// \f$ \int_0^T (r_\mathrm{des} - r) \mathrm{d}t \f$ [rad]
  double r_err_sum_;

  /// Accumulated heading error from controller
  /// \f$ \int_0^T (h_\mathrm{des} - h) \mathrm{d}t \f$ [rad * seconds]
  double h_err_sum_;

  /// Accumulated longitudinal velocity error from controller
  /// \f$ \int_0^T (u_\mathrm{des} - u) \mathrm{d}t \f$ [m]
  double u_err_sum_;

  // TODO check if u != w ever
  // TODO decide whether u_err_sum should be in msg

  /// Default constructor, zeroes all elements
  RoverState() : RoverState(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) {}

  /// Constructor
  /// \param x the x location of the rover
  /// \param y the y location of the rover
  /// \param h the heading of the rover
  /// \param u the longitudinal velocity of the rover
  /// \param v the latitudinal velocity of the rover
  /// \param r the yaw rate of the rover
  /// \param w the wheel speed of the rover
  RoverState(double x, double y, double h, double u, double v, double r,
             double w)
      : RoverState(x, y, h, u, v, r, w, 0.0, 0.0) {}

private:
  /// Constructor
  /// \param x the x location of the rover
  /// \param y the y location of the rover
  /// \param h the heading of the rover
  /// \param u the longitudinal velocity of the rover
  /// \param v the latitudinal velocity of the rover
  /// \param r the yaw rate of the rover
  /// \param w the wheel speed of the rover
  /// \param r_err_sum the accumulated r error during closed loop tracking
  /// \param h_err_sum the accumulated h error during closed loop tracking
  RoverState(double x, double y, double h, double u, double v, double r,
             double w, double r_err_sum, double h_err_sum)
      : x_{x}, y_{y}, h_{ToAngle(h)}, u_{u}, v_{v}, r_{r}, w_{w},
        r_err_sum_{r_err_sum}, h_err_sum_{h_err_sum},
        // TODO should u_err_sum always be zero?
        u_err_sum_{0.0} {}

public:
  /// Get the rover's heading
  /// \return the rover's heading in \f$ [-\pi, \pi] \f$
  double GetHeading() const;

  /// Set the rover's heading
  /// \param heading the heading to set the state to.
  void SetHeading(double heading);

  /// Generates a printable string of the state
  /// \return a string containing the x, y, h, u, v, and r state values
  std::string ToStringXYHUVR() const;

  /// Get the XYH components of the state
  /// \return the XYH components of the state
  PointXYH GetXYH() const;
};
} // namespace roahm
#endif
