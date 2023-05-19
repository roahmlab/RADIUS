#ifndef ROAHM_SIMPLE_UTIL_HPP_
#define ROAHM_SIMPLE_UTIL_HPP_

#include <algorithm> // for max, min
#include <cmath>     // for sqrt, cos, isfinite, sin
#include <complex>   // for arg, complex
#include <iostream>
#include <map>         // for map
#include <ostream>     // for operator<<, basic_ostream
#include <string>      // for operator<<, char_traits, string
#include <type_traits> // for declval, false_type, true_type
#include <utility>     // for swap

/// @file simple_util.hpp Contains many small utility functions

namespace roahm {

/// Converts a given input angle an angle in \f$[-\pi, \pi]\f$, returns 0 if
/// the input is not finite, for example if infinity or NaN's are provided.
/// \param theta the input angle
/// \return An equivalent angle to the input angle, restricted to
/// \f$[-\pi, \pi]\f$ iff the input argument is finite, otherwise returns 0
inline constexpr double ToAngle(double theta) {
  if (not std::isfinite(theta)) {
    return 0;
  }
  return std::arg(std::complex<double>{std::cos(theta), std::sin(theta)});
}

// Modified from https://stackoverflow.com/a/9154394
// https://stackoverflow.com/questions/257288/templated-check-for-the-
// existence-of-a-class-member-function
// (Xeo, Brad Larson)
namespace detail {
template <class> struct sfinae_true : std::true_type {};

template <class T, class A0>
static auto test_stream(int)
    -> sfinae_true<decltype(std::declval<T>().operator<<(std::declval<A0>()))>;
template <class, class A0> static auto test_stream(long) -> std::false_type;
} // namespace detail

template <class T, class Arg>
struct has_stream : decltype(detail::test_stream<T, Arg>(0)) {};

/// Clamps a value to be within a minimum and maximum. If the provided maximum
/// is less than the provided minimum, the minimum and maximum values will be
/// swapped prior to computation.
/// \param x the value to clamp
/// \param min_v the minimum value
/// \param max_v the maximum value
/// \return \f$ \begin{cases}
/// x & x \in [\textrm{min}, \textrm{max}] \\
/// \textrm{max} & x > \textrm{max} \\
/// \textrm{min} & x < \textrm{min}
/// \end{cases}
/// \f$
inline long Clamp(long x, long min_v, long max_v) {
  if (max_v < min_v) {
    std::swap(min_v, max_v);
  }
  return std::min(std::max(x, min_v), max_v);
}

/// Clamps a value to be within a minimum and maximum. The provided minimum and
/// maximum are expected to be in order.
/// \param x the value to clamp
/// \param min_v the minimum value
/// \param max_v the maximum value
/// \return \f$ \begin{cases}
/// x & x \in [\textrm{min}, \textrm{max}] \\
/// \textrm{max} & x > \textrm{max} \\
/// \textrm{min} & x < \textrm{min}
/// \end{cases}
/// \f$
inline double Clamp(double v, double min_v, double max_v) {
  // NaN returns false for all comparisons, so this
  // returns min_v if v is NaN
  if (v >= min_v and v <= max_v) {
    return v;
  } else if (v > max_v) {
    return max_v;
  } else if (v < min_v) {
    return min_v;
  }
  return min_v;
}

/// Clamps a value to be within a minimum and maximum. The provided minimum and
/// maximum are expected to be in order.
/// \param x the value to clamp
/// \param min_v the minimum value
/// \param max_v the maximum value
/// \param name the name of the parameter to use in warning
/// \return \f$ \begin{cases}
/// x & x \in [\textrm{min}, \textrm{max}] \\
/// \textrm{max} & x > \textrm{max} \\
/// \textrm{min} & x < \textrm{min}
/// \end{cases}
/// \f$
inline double ClampWithWarn(double x, double min_v, double max_v,
                            const std::string& name) {
  if (x > max_v or x < min_v) {
    std::cerr << name << " (" << x << ") not in [" << min_v << ", " << max_v
              << "]" << std::endl;
  }
  return Clamp(x, min_v, max_v);
}
/// Returns the square of the input
/// \param x the value to square
/// \return the square of the input, \f$ x^2 \f$
inline double Square(double x) { return x * x; }
/// Computes the L2 norm of the two input arguments as if they were a length-2
/// vector
/// \param x first element
/// \param y second element
/// \return \f$ \sqrt{x^2 + y^2} \f$
inline double Norm(double x, double y) { return std::sqrt(x * x + y * y); }

/// Checks whether the value is in the set \f$ [\textrm{min}, \textrm{max}) \f$,
/// assumes that the provided minimum and maximum are correctly ordered. If the
/// min or max are NaN, this will always return false.
/// \param x the value to check
/// \param min_v the minimum of the set, inclusive
/// \param max_v the maximum of the set, exclusive
/// \return true iff \f$ x \in [\textrm{min}, \textrm{max}) \f$, false
/// otherwise, or if min or max are NaN values
inline bool InClosedOpenInterval(double x, double min_v, double max_v) {
  // Returns true iff val in [min_v, max_v)
  // NaN min or max always return false
  return (min_v <= x) and (x < max_v);
}

/// Checks whether the value is in the set \f$ [\textrm{min}, \textrm{max}] \f$,
/// assumes that the provided minimum and maximum are correctly ordered. If the
/// min or max are NaN, this will always return false.
/// \param x the value to check
/// \param min_v the minimum of the set, inclusive
/// \param max_v the maximum of the set, inclusive
/// \return true iff \f$ x \in [\textrm{min}, \textrm{max}] \f$, false
/// otherwise, or if min or max are NaN values
inline bool InClosedInterval(double x, double min_v, double max_v) {
  // Returns true iff val in [min_v, max_v]
  // NaN min or max always return false
  return (min_v <= x) and (x <= max_v);
}

/// Computes the shortest distance between two angles.
/// \param theta1 the angle to "end" at [rad]
/// \param theta2 the angle to "start" at [rad]
/// \return \f$ \theta_1 - \theta_2 \f$ radians, in the range
/// \f$ [-\pi, \pi] \f$, returns 0 if either of the arguments are NaN values
inline constexpr double AngleDiff(double theta1, double theta2) { // JL CHECKED
  return ToAngle(theta1 - theta2);
}

/// Represents an interval \f$ [\textrm{min}, \textrm{max}] \f$
class Interval {
  /// The maximum value of the interval
  double min_;
  /// The minimum value of the interval
  double max_;

public:
  /// Default constructor
  Interval() : min_{0.0}, max_{0.0} {}

  /// Constructor. If the provided maximum is less than the provided minimum,
  /// the bounds will be swapped.
  /// \param min_val The minimum value of the interval
  /// \param max_val The maximum value of the interval
  Interval(double min_val, double max_val) : min_{0}, max_{0} {
    SetMin(min_val);
    SetMax(max_val);
  }

  /// Computes the distance, \f$ d \f$ from the interval to some number
  /// \f$ x \f$, where the distance is defined as:
  /// \f$ d = \begin{cases}
  /// 0 & x \in [\textrm{min}, \textrm{max}] \\
  /// x - \textrm{max} & x > \textrm{max} \\
  /// \textrm{min} - x & x < \textrm{min}
  /// \end{cases} \f$
  /// \param x the number to compute the distance to
  /// \return the distance from the interval to the input, as defined above
  double DistanceTo(double x) const;

  /// Gets the minimum value of the interval
  /// \return the minimum value of the interval
  inline double Min() const { return min_; }

  /// Gets the maximum value of the interval
  /// \return the maximum value of the interval
  inline double Max() const { return max_; }

  /// Sets the minimum value of the interval. If the current maximum is less
  /// than the new minimum, the maximum and minimum will be set to the new value
  /// \param min_val the new minimum value of the interval
  inline void SetMin(double min_val) {
    min_ = min_val;
    if (max_ < min_) {
      max_ = min_;
    }
  }

  /// Sets the maximum value of the interval. If the new maximum is less than
  /// the current minimum, the maximum and minimum will be set to the new value
  /// \param max_val the value to set the maximum end of the interval to.
  inline void SetMax(double max_val) {
    max_ = max_val;
    if (max_ < min_) {
      min_ = max_;
    }
  }

  // Gets the interval width \f$ | \textrm{max} - \textrm{min} | \f$
  // \return the interval width \f$ | \textrm{max} - \textrm{min} | \f$
  inline double Width() const { return std::abs(max_ - min_); }

  inline double Midpoint() const { return min_ + ((max_ - min_) * 0.5); }

  inline bool operator==(const Interval& other) const {
    return (min_ == other.min_) and (max_ == other.max_);
  }
};

/// Returns a comma-separated string "arg0, arg1, ..., argN" using the << on all
/// provided arguments
/// \tparam T type of the first argument
/// \tparam U type list
/// \param arg0 the first element of the list to print
/// \param args remaining arguments to print
/// \return a comma-separated string "arg0, arg1, ..., argN"
template <typename T, typename... U>
std::string VariadicListToStr(T arg0, U... args) {
  std::stringstream s;
  s << arg0;
  ((s << ", " << args), ...);
  return s.str();
}

} // namespace roahm
#endif
