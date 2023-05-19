#include "risk_constraint.hpp"

#include <chrono>
#include <cmath>
#include <stdexcept>

#include "manu_type.hpp"
#include "timing_util.hpp"
#include "unit_conversion.hpp"
#include <fmt/format.h>

namespace roahm::risk_constraint {
namespace {

double IntPow(const double val, int exp_val) {
  double ret = 1.0;
  const auto abs_exp_val = std::abs(exp_val);
  for (int i = 0; i < abs_exp_val; ++i) {
    ret *= val;
  }

  if (exp_val < 0) {
    ret = 1.0 / ret;
  }
  return ret;
}

double RiskThresholdAt(const double u_meters_per_second) {
  // Polynomial fit to data in NHTSA paper (TODO link?)
  const double u_mph = MetersPerSecondToMilesPerHour(u_meters_per_second);
  return 0.01 * (1.53e-6 * IntPow(u_mph, 4) - 1.84e-4 * IntPow(u_mph, 3) +
                 5.49e-3 * IntPow(u_mph, 2) + 0.011 * (u_mph) + 2.23e-3);
}
double RiskThresholdGradientAt(const double u_meters_per_second) {
  const double u_mph = MetersPerSecondToMilesPerHour(u_meters_per_second);
  return 0.01 * MetersPerSecondToMilesPerHour(1.53e-6 * IntPow(u_mph, 3) * 4 -
                                              1.84e-4 * IntPow(u_mph, 2) * 3 +
                                              5.49e-3 * (u_mph)*2 + 0.011);
}
} // namespace

RiskConstraint::RiskThresholdAndGradient
ComputeRiskThresholdAndGradient(const double u_meters_per_sec,
                                const bool is_spd, const double param_val) {
  const double tracking_error = 0.0;
  const double u0_min = u_meters_per_sec - tracking_error;
  const double u0_max = u_meters_per_sec + tracking_error;
  const double pmt = param_val - tracking_error;
  const double ppt = param_val + tracking_error;
  RiskConstraint::RiskThresholdAndGradient risk_thresh_and_grad_out{0.0, 0.0};
  // TODO better name
  constexpr double kCrit = MilesPerHourToMetersPerSecond(59.574);
  const bool is_constant_spd = (not is_spd) or (param_val != u_meters_per_sec);
  if (not is_constant_spd) {
    if (pmt < u0_min) {
      if (u0_max < kCrit) {
        // L597
        const auto th1 = RiskThresholdAt(pmt);
        const auto th2 = RiskThresholdAt(u0_max);
        if (th1 < th2) {
          risk_thresh_and_grad_out = {th1, RiskThresholdGradientAt(pmt)};
        } else {
          risk_thresh_and_grad_out = {th2, 0.0};
        }
      } else if (pmt < kCrit) {
        // L609
        risk_thresh_and_grad_out = {RiskThresholdAt(kCrit), 0.0};
      } else {
        // L613
        risk_thresh_and_grad_out = {RiskThresholdAt(pmt),
                                    RiskThresholdGradientAt(pmt)};
      }
    } else /* if (ppt > u0_max) */ {
      // L618
      if (ppt < kCrit) {
        // L619
        // NOTE TODO differs from MATLAB since MATLAB presumably
        // has bug (L619 as of 2d64931)
        const auto th1 = RiskThresholdAt(ppt);
        const auto th2 = RiskThresholdGradientAt(u0_min);
        if (th1 < th2) {
          // L622
          risk_thresh_and_grad_out = {th1, RiskThresholdGradientAt(ppt)};
        } else {
          // L626
          risk_thresh_and_grad_out = {th2, 0.0};
        }
      } else if (u0_min < kCrit) {
        // L631
        risk_thresh_and_grad_out = {RiskThresholdAt(kCrit), 0.0};
      } else {
        // L635
        risk_thresh_and_grad_out = {RiskThresholdAt(u0_max), 0.0};
      }
    }
  } else {
    // Constant speed
    // NOTE (Maybe TODO):
    // the comparisons (p +/- tracking) < (or >) (u0 +/- tracking) are
    // equivalent to p < u0 and p > u0, which is why even in the case of a speed
    // change, this is constant speed (and why there is no else in the if speed
    // section)

    // L640
    if (u0_max < kCrit) {
      // L641
      const auto th1 = RiskThresholdAt(u0_min);
      const auto th2 = RiskThresholdAt(u0_max);
      risk_thresh_and_grad_out = {std::min(th1, th2), 0.0};
    } else if (u0_min < kCrit) {
      // L648
      risk_thresh_and_grad_out = {RiskThresholdAt(kCrit), 0.0};
    } else {
      // L652
      risk_thresh_and_grad_out = {RiskThresholdAt(u0_min), 0.0};
    }
  }
  // fmt::print("Risk Threshold: {}\n",
  // risk_thresh_and_grad_out.risk_threshold_);
  return risk_thresh_and_grad_out;
}
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
  fmt::print("Pre Slicing took: {}\n", GetDeltaS(t1, t0));
}
RiskConstraintValues RiskConstraint::EvaluateAt(
    double p_in,
    const RiskThresholdAndGradient risk_thresh_and_gradient) noexcept(false) {
  // fmt::print("Risk Constraint Evaluation\n");
  //  The output risk and its gradient and hessian
  RiskConstraintValues out{};

  // The parameter lambda (sliced)
  const double p_val = (p_in - cg_p_.at(0)) / cg_p_.at(1);

  // Sum the risk (and its gradients and hessians) for all obstacles
  for (auto prob_inputs : all_prob_inputs_) {
    prob_inputs.p_ = p_val;
    // const auto t0 = Tick();
    const auto prob_outputs = ::roahm::prob_integ::ProbIntegration(prob_inputs);
    // const auto t1 = Tick();
    // fmt::print(" Prob Integ []: {}\n", prob_outputs.constraint_val_);
    //  const std::string out_str = "ProbIntegration took " +
    //                              std::to_string(GetDeltaS(t1, t0) * 1.0e3) +
    //                              "ms\n";
    //  std::cout << out_str << std::flush;
    out.risk_ += prob_outputs.constraint_val_;
    out.d_risk_ += prob_outputs.d_constraint_val_;
    out.d2_risk_ += prob_outputs.d2_constraint_val_;
    // fmt::print("Prob Outputs Cons Val:   {}\n",
    // prob_outputs.constraint_val_); fmt::print("Prob Outputs DCons Val: {}\n",
    // prob_outputs.d_constraint_val_); fmt::print("Prob Outputs D2Cons Val:
    // {}\n",
    //            prob_outputs.d2_constraint_val_);
  }

  // Since our constraint is in the form Risk < epsilon,
  // we can instead use (Risk - epsilon) < 0 and keep the
  // constraint upper bound constant at zero in the optimization problem
  // fmt::print(" Risk Threshold: {}\n",
  // risk_thresh_and_gradient.risk_threshold_);
  // fmt::print(" Risk Thresh:   {}\n",
  // risk_thresh_and_gradient.risk_threshold_); fmt::print(" Risk Thresh G:
  // {}\n",
  //           risk_thresh_and_gradient.risk_threshold_gradient_);
  out.risk_ -= risk_thresh_and_gradient.risk_threshold_;
  out.d_risk_ -= risk_thresh_and_gradient.risk_threshold_gradient_;

  // fmt::print(" Out Risk Thresh:  {}\n",
  //            risk_thresh_and_gradient.risk_threshold_);
  // fmt::print(" Out Risk Cons:    {}\n", out.risk_);
  // fmt::print(" Out D Risk Cons:  {}\n", out.d_risk_);
  // fmt::print(" Out D2 Risk Cons: {}\n", out.d2_risk_);

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
RiskConstraintValues RiskConstraint::EvaluateAt(double p_in) noexcept(false) {
  const auto risk_threshold_and_gradient = ComputeRiskThresholdAndGradient(
      u0_meters_per_second_, is_speed_change_manu_, p_in);
  return EvaluateAt(p_in, risk_threshold_and_gradient);
}
} // namespace roahm::risk_constraint