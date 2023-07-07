#include "risk_constraint.hpp"

#include <chrono>
#include <cmath>
#include <stdexcept>

#include "manu_type.hpp"
#include "timing_util.hpp"
#include "unit_conversion.hpp"
#include <fmt/format.h>

namespace roahm::risk_constraint {

RiskConstraint::RiskConstraint(
    const Vehrs& vehrs, const std::array<double, 3>& u0v0r0_slice_beta,
    const std::vector<std::vector<double>>& all_mu_sigmas,
    const double u0_meters_per_second)
    : cuda_info_{vehrs.cuda_info_}, all_mu_sigmas_{all_mu_sigmas},
      all_pre_slice_outputs_{}, cg_p_{cuda_info_.cg_p_.at(0),
                                      cuda_info_.cg_p_.at(1)},
      all_prob_inputs_{}, u0_meters_per_second_{u0_meters_per_second},
      is_speed_change_manu_{IsSpd(vehrs.GetManuType())} {
  const auto t0 = Tick();
  for (const auto& mu_sigma_individ : all_mu_sigmas_) {
    all_pre_slice_outputs_.push_back(::roahm::pre_slice::PreSlice(
        cuda_info_, u0v0r0_slice_beta, mu_sigma_individ.data()));
    const auto& pre_slice_outputs = all_pre_slice_outputs_.back();
    all_prob_inputs_.push_back(::roahm::ProbIntegrationInputs{
        cuda_info_.num_zono_, cuda_info_.grid_x0_.data(),
        cuda_info_.grid_y0_.data(), cuda_info_.grid_dx_.data(),
        cuda_info_.grid_dy_.data(), cuda_info_.block_inzono_list_.data(),
        cuda_info_.rot_angle_.data(), mu_sigma_individ.data(),
        cuda_info_.cg_p_.data(), cuda_info_.g_p_x_.data(),
        cuda_info_.grid_size_, 0.0, pre_slice_outputs.H1_.get(),
        pre_slice_outputs.H2_.get(), pre_slice_outputs.H4_.get()});
  }
  const auto t1 = Tick();
}
RiskConstraintValues RiskConstraint::EvaluateAt(
    double p_in,
    const RiskThresholdAndGradient risk_thresh_and_gradient) noexcept(false) {
  //  The output risk and its gradient and hessian
  RiskConstraintValues out{};

  // The parameter lambda (sliced)
  const double p_val = (p_in - cg_p_.at(0)) / cg_p_.at(1);

  // Sum the risk (and its gradients and hessians) for all obstacles
  for (auto prob_inputs : all_prob_inputs_) {
    prob_inputs.p_ = p_val;
    const auto prob_outputs = ::roahm::prob_integ::ProbIntegration(prob_inputs);
    out.risk_ += prob_outputs.constraint_val_;
    out.d_risk_ += prob_outputs.d_constraint_val_;
    out.d2_risk_ += prob_outputs.d2_constraint_val_;
  }

  // Since our constraint is in the form Risk < epsilon,
  // we can instead use (Risk - epsilon) < 0 and keep the
  // constraint upper bound constant at zero in the optimization problem
  out.risk_ -= risk_thresh_and_gradient.risk_threshold_;
  out.d_risk_ -= risk_thresh_and_gradient.risk_threshold_gradient_;

  // Check if any of the CUDA outputs or our risk threshold are NaN
  if (std::isnan(risk_thresh_and_gradient.risk_threshold_)) {
    throw std::runtime_error("Risk threshold is nan");
  }

  if (std::isnan(risk_thresh_and_gradient.risk_threshold_gradient_)) {
    throw std::runtime_error("Risk threshold gradient is nan");
  }

  if (std::isnan(out.risk_)) {
    throw std::runtime_error("Risk constraint is nan");
  }
  if (std::isnan(out.d_risk_)) {
    throw std::runtime_error("Risk constraint gradient is nan");
  }
  if (std::isnan(out.d2_risk_)) {
    throw std::runtime_error("Risk constraint hessian is nan");
  }

  return out;
}

RiskConstraintValues
RiskConstraint::EvaluateAt(double p_in,
                           const double fixed_risk_threshold) noexcept(false) {
  return EvaluateAt(p_in, {fixed_risk_threshold, 0.0});
}
} // namespace roahm::risk_constraint
