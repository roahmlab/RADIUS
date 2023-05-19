#ifndef ROAHM_RISK_CONSTRAINT_HPP_
#define ROAHM_RISK_CONSTRAINT_HPP_
#include <array>
#include <vector>

#include "cuda_info.hpp"
#include "pre_slice.hpp"
#include "prob_integ.hpp"
#include "prob_integ_io.hpp"
#include "vehrs.hpp"
namespace roahm::risk_constraint {

double ComputeRiskThreshold(const double u_meters_per_sec);

struct RiskConstraintValues {
  double risk_;
  double d_risk_;
  double d2_risk_;
};

struct RiskConstraint {
public:
  struct RiskThresholdAndGradient {
    double risk_threshold_;
    double risk_threshold_gradient_;
  };

private:
  // TODO maybe don't copy mu and cuda info eventually
  const ::roahm::CudaInfo cuda_info_;
  std::vector<std::vector<double>> all_mu_sigmas_;
  std::vector<::roahm::pre_slice::PreSliceOutputs> all_pre_slice_outputs_;
  const std::array<double, 2> cg_p_;
  std::vector<::roahm::ProbIntegrationInputs> all_prob_inputs_;
  const double u0_meters_per_second_;
  const bool is_speed_change_manu_;
public:
  RiskConstraintValues EvaluateAt(
      double p_in,
      const RiskThresholdAndGradient risk_thresh_and_gradient) noexcept(false);

  RiskConstraint(const Vehrs& vehrs,
                 const std::array<double, 3>& u0v0r0_slice_beta,
                 const std::vector<std::vector<double>>& all_mu_sigmas,
                 const double u0_meters_per_second);

  RiskConstraintValues EvaluateAt(double p_in) noexcept(false);

  /// @brief Evaluate the risk constraint, and its derivatives at a given point,
  /// with some maximum threshold
  /// @param p_in the trajectory parameter to evaluate the constraint at
  /// @param fixed_risk_threshold the threshold to use as an upper bound on the
  /// constraint
  /// @return the risk constraint minus the risk threshold, and the gradients of
  /// the risk
  RiskConstraintValues
  EvaluateAt(double p_in, const double fixed_risk_threshold) noexcept(false);

  [[nodiscard]] int GetNumMuSigmas() const noexcept {
    return all_mu_sigmas_.size();
  }
  [[nodiscard]] int GetMuSigmaWidth() const noexcept {
    return all_mu_sigmas_.empty() ? 0 : (all_mu_sigmas_.at(0).size());
  }
};

} // namespace roahm::risk_constraint
#endif
