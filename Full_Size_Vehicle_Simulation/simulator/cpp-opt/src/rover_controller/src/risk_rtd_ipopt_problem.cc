#include "risk_rtd_ipopt_problem.hpp"

#include <chrono>
#include <cmath>
#include <fmt/format.h>
#include <iostream>
#include <string>

#include "IpAlgTypes.hpp"
#include "IpTNLP.hpp"
#include "fl_zono_obs_set.hpp"
#include "gencon.hpp"
#include "ipopt_string_utils.hpp"
#include "manu_type.hpp"
#include "risk_constraint.hpp"
#include "timing_util.hpp"
#include "waypoint.hpp"

namespace roahm::risk_rtd_ipopt_problem {

namespace {
using Number = Ipopt::Number;
using Index = Ipopt::Index;

} // namespace

std::vector<Ipopt::SmartPtr<Ipopt::IpoptApplication>>
GetRiskIpoptApps(int num_apps) {
  // param: number of ipopt you want to construct, equals to number of threads
  // returns vector of Ipopt solvers and init their params

  std::vector<Ipopt::SmartPtr<Ipopt::IpoptApplication>> app_vec;
  app_vec.reserve(num_apps);
  for (int i = 0; i < num_apps; ++i) {
    app_vec.push_back(IpoptApplicationFactory());
  }
  for (auto& app : app_vec) {
    // Refer to https://www.gams.com/latest/docs/S_IPOPT.html
    auto opts = app->Options();
    opts->SetNumericValue("tol", 1.0e-6);
    opts->SetStringValue("hessian_constant", "no");
    opts->SetStringValue("linear_solver", "ma57");
    opts->SetNumericValue("ma57_pre_alloc", 5.0);
    opts->SetIntegerValue("print_level", 1);
    opts->SetNumericValue("max_wall_time", 1.0 / 5.0);
    opts->SetIntegerValue("max_iter", 25);
    opts->SetStringValue("jacobian_approximation", "finite-difference-values");
  }
  return app_vec;
}

RiskRtdIpoptProblem::RiskRtdIpoptProblem(
    const Vehrs& vehrs, const std::array<double, 3>& u0v0r0_slice_beta,
    const OptimizationInputs& opt_inputs)
    : // Risk constraints
      risk_constraint_enabled_{opt_inputs.HasRiskConstraints()},
      num_risk_constraints_{opt_inputs.NumRiskConstraints()},
      risk_constraint_{vehrs, u0v0r0_slice_beta, opt_inputs.mu_sigmas_,
                       opt_inputs.GetU0MetersPerSecond()},
      last_risk_constraint_values_{},
      risk_threshold_{opt_inputs.risk_threshold_},
      // Decision variable bounds
      trajectory_param_min_{vehrs.GetTrajParamMin()},
      trajectory_param_max_{vehrs.GetTrajParamMax()},
      // Classical FL Zono constraints
      classical_constraints_{::roahm::fl_zono_constraint::GetFlZonoConstraint(
          vehrs, opt_inputs.fl_zono_obs_)},
      num_classical_constraints_{classical_constraints_.NumConstraints()},
      classical_constraint_enabled_{opt_inputs.HasClassicalConstraints()},
      // Cached values
      last_trajectory_param_{std::nan("1")}, last_cost_values_{},
      // Cost function
      waypoint_cost_{WaypointCostFromVehrsAndWp(
          vehrs, opt_inputs.waypoint_heuristic_adjusted_, opt_inputs.use_left_cost_function_)},
      timings_{}, final_param_{std::nullopt}, final_cost_{std::nullopt},
      final_location_{std::nullopt},
      waypoint_local_mirror_accounted_{
          opt_inputs.waypoint_local_mirror_accounted_},
      waypoint_heuristic_adjusted_{opt_inputs.waypoint_heuristic_adjusted_},
      found_feasible_internal_{false} {
  // Classical FL Zono Constraints

  if ((not risk_constraint_enabled_) and (num_risk_constraints_ != 0)) {
    throw std::runtime_error("Risk constraints disabled but num > 0!");
  }
  if (risk_constraint_enabled_ and (num_risk_constraints_ != 1)) {
    throw std::runtime_error("Risk constraints enabled but num != 1");
  }

  if (not opt_inputs.use_fixed_risk_threshold_) {
    throw std::runtime_error("Currently only fixed risk threshold is enabled");
  }
}

void RiskRtdIpoptProblem::Recompute(const double trajectory_param) {
  if (last_trajectory_param_ != trajectory_param) {
    if (risk_constraint_enabled_) {
      last_trajectory_param_ = trajectory_param;
      const auto t0 = Tick();
      const auto risk_threshold_and_gradient{::roahm::risk_constraint::RiskConstraint::RiskThresholdAndGradient{risk_threshold_, 0.0}};
      last_risk_constraint_values_ = risk_constraint_.EvaluateAt(trajectory_param, risk_threshold_and_gradient);
      const auto t1 = Tick();
      timings_.risk_evals_seconds_ += GetDeltaS(t1, t0);
      ++timings_.num_risk_evals_;
      last_cost_values_ = waypoint_cost_.EvaluateAt(trajectory_param);
    }
  }
}

bool RiskRtdIpoptProblem::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                       Index& nnz_h_lag,
                                       IndexStyleEnum& index_style) {
  // Number of decision variables
  n = 1;
  // Number of constraints
  m = num_classical_constraints_ + num_risk_constraints_;

  // Nonzero jacobian entries
  nnz_jac_g = n * m;

  // Nonzero hessian of lagrangian
  nnz_h_lag = 1;
  index_style = TNLP::C_STYLE;
  return true;
}
bool RiskRtdIpoptProblem::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                          Index m, Number* g_l, Number* g_u) {
  x_l[0] = trajectory_param_min_;
  x_u[0] = trajectory_param_max_;

  if (m != num_classical_constraints_ + num_risk_constraints_) {
    throw std::runtime_error("m != num classical + num risk");
  }

  for (Index r = 0; r < num_classical_constraints_; ++r) {
    g_l[r] = -1.0e+19;
    g_u[r] = 0.0;
  }

  for (Index r = num_classical_constraints_;
       r < num_classical_constraints_ + num_risk_constraints_; ++r) {
    g_l[r] = -1.0e+19;
    g_u[r] = 0.0;
  }

  return true;
}
bool RiskRtdIpoptProblem::get_starting_point(Index n, bool init_x, Number* x,
                                             bool init_z, Number* z_L,
                                             Number* z_U, Index m,
                                             bool init_lambda, Number* lambda) {
  if (not init_x) {
    throw std::runtime_error("init_x must be true");
  }
  if (init_z) {
    throw std::runtime_error("init_z must be false");
  }
  if (init_lambda) {
    throw std::runtime_error("init_lambda must be false");
  }

  x[0] = (trajectory_param_min_ + trajectory_param_max_) / 2.0;
  return true;
}

bool RiskRtdIpoptProblem::eval_f(Index n, const Number* x, bool new_x,
                                 Number& obj_value) {
  Recompute(x[0]);
  obj_value = last_cost_values_.cost_;
  return true;
}

bool RiskRtdIpoptProblem::eval_grad_f(Index n, const Number* x, bool new_x,
                                      Number* grad_f) {
  Recompute(x[0]);
  grad_f[0] = last_cost_values_.cost_gradient_wrt_param_;
  return true;
}

bool RiskRtdIpoptProblem::eval_g(Index n, const Number* x, bool new_x, Index m,
                                 Number* g) {
  Recompute(x[0]);
  bool found_cons_viol = false;
  if (classical_constraint_enabled_) {
    classical_constraints_.GetConstraintEvaluations(x[0], g,
                                                    num_classical_constraints_);
    for (int i = 0; i < num_classical_constraints_; ++i) {
      if (g[i] >= 0.0) {
        found_cons_viol = true;
        break;
      }
    }
  }
  if (risk_constraint_enabled_) {
    g[num_classical_constraints_] = last_risk_constraint_values_.risk_;
    if (g[num_classical_constraints_] > 0) {
      found_cons_viol = true;
    }
  }
  if (not found_cons_viol) {
    last_cost_values_ = waypoint_cost_.EvaluateAt(x[0]);
    if ((not final_cost_.has_value()) or
        (last_cost_values_.cost_ < final_cost_.value())) {
      found_feasible_internal_ = true;
      final_param_ = x[0];
      final_cost_ = last_cost_values_.cost_;
      final_location_ = waypoint_cost_.GetPointAt(x[0]);
    }
  }
  return true;
}

bool RiskRtdIpoptProblem::eval_jac_g(Index n, const Number* x, bool new_x,
                                     Index m, Index nele_jac, Index* iRow,
                                     Index* jCol, Number* values) {
  if (values == NULL) {
    // Return structure of the jacobian of the constraints
    // element at i,j: grad_{x_j} g_{i}(x)
    for (Index r = 0; r < m; ++r) {
      for (Index c = 0; c < n; ++c) {
        Index idx = (r * n) + c;
        iRow[idx] = r;
        jCol[idx] = c;
      }
    }
  } else {
    Recompute(x[0]);
    if (classical_constraint_enabled_) {
      classical_constraints_.GetGradient(x[0], values,
                                         num_classical_constraints_);
    }
    if (risk_constraint_enabled_) {
      values[num_classical_constraints_] = last_risk_constraint_values_.d_risk_;
    }
  }
  return true;
}

bool RiskRtdIpoptProblem::eval_h(Index n, const Number* x, bool new_x,
                                 Number obj_factor, Index m,
                                 const Number* lambda, bool new_lambda,
                                 Index nele_hess, Index* iRow, Index* jCol,
                                 Number* values) {
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // element at 1,1: grad^2_{x1,x1} L(x,lambda)
    iRow[0] = 0;
    jCol[0] = 0;
  } else {
    Recompute(x[0]);
    // return the values
    // element at 1,1: grad^2_{x1,x1} L(x,lambda)
    values[0] = obj_factor * last_cost_values_.cost_hessian_wrt_param_;
    if (risk_constraint_enabled_) {
      values[0] += lambda[0] * last_risk_constraint_values_.d2_risk_;
    }
  }
  return true;
}

void RiskRtdIpoptProblem::finalize_solution(
    Ipopt::SolverReturn status, Index n, const Number* x, const Number* z_L,
    const Number* z_U, Index m, const Number* g, const Number* lambda,
    Number obj_value, const Ipopt::IpoptData* ip_data,
    Ipopt::IpoptCalculatedQuantities* ip_cq) {
  final_return_status_ = status;
  if (status == Ipopt::SolverReturn::SUCCESS) {
    found_feasible_internal_ = true;
    final_param_ = x[0];
    final_cost_ = obj_value;
    final_location_ = waypoint_cost_.GetPointAt(x[0]);
  }
}

bool RiskRtdIpoptProblem::FoundFeasible() const {
  return found_feasible_internal_;
}

std::optional<Ipopt::SolverReturn>
RiskRtdIpoptProblem::GetSolverReturnStatus() const {
  return final_return_status_;
}

std::optional<double> RiskRtdIpoptProblem::GetFeasibleCost() const {
  return final_cost_;
}

std::optional<double> RiskRtdIpoptProblem::GetFeasibleParam() const {
  return final_param_;
}
std::optional<PointXYH> RiskRtdIpoptProblem::GetFeasibleLocation() const {
  return final_location_;
}

WaypointLocalMirrorTakenIntoAccount
RiskRtdIpoptProblem::GetWaypointHeuristicAdjusted() const {
  return waypoint_heuristic_adjusted_;
}

WaypointLocalMirrorTakenIntoAccount
RiskRtdIpoptProblem::GetWaypointLocalMirrorAccounted() const {
  return waypoint_local_mirror_accounted_;
}

} // namespace roahm::risk_rtd_ipopt_problem
