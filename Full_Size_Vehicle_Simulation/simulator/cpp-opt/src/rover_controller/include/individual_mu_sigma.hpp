#ifndef ROAHM_INDIVIDUAL_MU_SIGMA_HPP_
#define ROAHM_INDIVIDUAL_MU_SIGMA_HPP_

#include "point_xyh.hpp"

namespace roahm {

/// @brief Contains one set of mu and sigma values representing an obstacle
/// probability density function (PDF), along with a helper function to
/// convert it from its current frame into another frame of reference.
struct IndividualMuSigma {
  /// @brief The mean x location
  const double mu_x_;

  /// @brief The mean y location
  const double mu_y_;

  /// @brief The covariance matrix's x-x covariance. (NOT standard deviation)
  const double sigma_1_;

  /// @brief The covariance matrix's x-y covariance. (NOT standard deviation)
  const double sigma_2_;

  /// @brief The covariance matrix's y-y covariance. (NOT standard deviation)

  const double sigma_4_;

  /// @brief Compute the same PDF but relative to a different frame of reference
  /// @param frame the new frame of reference
  /// @param mirror whether or not the PDF is going to be mirrored about the
  /// x-axis
  /// @return the new PDF relative to @p frame and, if @p mirror is true,
  /// mirrored about the x-axis
  IndividualMuSigma RelativeTo(const ::roahm::PointXYH& frame,
                               const bool mirror) const;
};
} // namespace roahm

#endif // ROAHM_INDIVIDUAL_MU_SIGMA_HPP_