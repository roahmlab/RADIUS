#include "monte_carlo_ipopt.hpp"
#include "comparison_methods.hpp"
#include "ipopt_string_utils.hpp"
#include "rover_state.hpp"
#include "timing_util.hpp"
#include <random>
#include <stdexcept>

namespace roahm::monte_carlo_ipopt {

MonteCarloIpoptProblem::MonteCarloIpoptProblem(
    const Vehrs& vehrs, const std::array<double, 3>& u0v0r0_slice_beta,
    const OptimizationInputs& opt_inputs,
    // New
    const std::size_t seed_val, const Vehrs& unsliced_vehrs,
    const RoverState& z0, const risk_comparisons::EnvFootprints footprints,
    const int num_trajectory_samples)
    : // Risk constraints
      risk_constraint_enabled_{opt_inputs.HasRiskConstraints()},
      num_risk_constraints_{opt_inputs.NumRiskConstraints()},
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
          vehrs, opt_inputs.waypoint_heuristic_adjusted_)},
      timings_{}, gen_{seed_val}, final_param_{std::nullopt},
      final_cost_{std::nullopt}, final_location_{std::nullopt},
      waypoint_local_mirror_accounted_{
          opt_inputs.waypoint_local_mirror_accounted_},
      waypoint_heuristic_adjusted_{opt_inputs.waypoint_heuristic_adjusted_},
      found_feasible_internal_{false}, unsliced_vehrs_{unsliced_vehrs}, z0_{z0},
      footprints_{footprints}, all_mu_sigmas_{opt_inputs.mu_sigmas_},
      num_trajectory_samples_{num_trajectory_samples},
      last_risk_constraint_values_{0.0} {
  // Classical FL Zono Constraints
  // classical_constraint_enabled_ = opt_inputs.HasClassicalConstraints();
  // if (classical_constraint_enabled_) {
  //  num_classical_constraints_ =
  //  classical_constraints_.value().zono_startpoints_.size();
  //}

  if ((not risk_constraint_enabled_) and (num_risk_constraints_ != 0)) {
    throw std::runtime_error("Risk constraints disabled but num > 0!");
  }
  if (risk_constraint_enabled_ and (num_risk_constraints_ != 1)) {
    throw std::runtime_error("Risk constraints enabled but num != 1");
  }

  fmt::print("[RISK IPOPT] Classical Constraints Enabled: {}\n",
             classical_constraint_enabled_);
  fmt::print("[RISK IPOPT] Num Classical Constraints:     {}\n",
             num_classical_constraints_);
  fmt::print("[RISK IPOPT] Risk Constraints Enabled:      {}\n",
             risk_constraint_enabled_);
  fmt::print("[RISK IPOPT] Num Risk Constraints:          {}\n",
             num_risk_constraints_);
}

void MonteCarloIpoptProblem::Recompute(const double trajectory_param) {
  if (last_trajectory_param_ != trajectory_param) {
    if (risk_constraint_enabled_) {
      last_trajectory_param_ = trajectory_param;
      const auto t0 = Tick();
      const double risk_threshold{0.05};
      const Vehrs param_sliced_frs{unsliced_vehrs_.SliceAtParam(
          z0_.u_, z0_.v_, z0_.r_, trajectory_param)};
      const auto opt_e_con{risk_comparisons::OptEConstraint(
          param_sliced_frs, all_mu_sigmas_, footprints_,
          num_trajectory_samples_, risk_threshold, gen_)};
      last_risk_constraint_values_ =
          opt_e_con.risk_ - opt_e_con.risk_threshold_;
      const auto t1 = Tick();
      timings_.risk_evals_seconds_ += GetDeltaS(t1, t0);
      ++timings_.num_risk_evals_;
      last_cost_values_ = waypoint_cost_.EvaluateAt(trajectory_param);
      /*
      std::cout << "Param:    " << trajectory_param << " of [" <<
      trajectory_param_min_ << ", " << trajectory_param_max_ << "]" <<
      std::endl; std::cout << " Risk:    " << last_risk_constraint_values_.risk_
      << std::endl; std::cout << " D Risk:  " <<
      last_risk_constraint_values_.d_risk_ << std::endl; std::cout << " D2 Risk:
      " << last_risk_constraint_values_.d2_risk_ << std::endl; std::cout << "
      Cost:    " << last_cost_values_.cost_ << std::endl; std::cout << " D Cost:
      " << last_cost_values_.cost_gradient_wrt_param_ << std::endl; std::cout <<
      " D2 Cost: " << last_cost_values_.cost_hessian_wrt_param_<< std::endl;
      */
    }
  }
}

bool MonteCarloIpoptProblem::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
                                          Ipopt::Index& nnz_jac_g,
                                          Ipopt::Index& nnz_h_lag,
                                          IndexStyleEnum& index_style) {
  // Ipopt::Number of decision variables
  n = 1;
  // Ipopt::Number of constraints
  m = num_classical_constraints_ + num_risk_constraints_;

  // Nonzero jacobian entries
  nnz_jac_g = n * m;

  // Nonzero hessian of lagrangian
  nnz_h_lag = 1;
  index_style = TNLP::C_STYLE;
  return true;
}
bool MonteCarloIpoptProblem::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l,
                                             Ipopt::Number* x_u, Ipopt::Index m,
                                             Ipopt::Number* g_l,
                                             Ipopt::Number* g_u) {
  // fmt::print("Trajectory Range [{}, {}]\n", trajectory_param_min_,
  //            trajectory_param_max_);
  x_l[0] = trajectory_param_min_;
  x_u[0] = trajectory_param_max_;

  if (m != num_classical_constraints_ + num_risk_constraints_) {
    throw std::runtime_error("m != num classical + num risk");
  }

  // TODO
  for (Ipopt::Index r = 0; r < num_classical_constraints_; ++r) {
    g_l[r] = -1.0e+19;
    g_u[r] = 0.0;
  }

  // TODO
  for (Ipopt::Index r = num_classical_constraints_;
       r < num_classical_constraints_ + num_risk_constraints_; ++r) {
    // TODO should this just be -inf?
    // Set this to extremely low number, since
    g_l[r] = -1.0e+19;
    g_u[r] = 0.0;
  }

  return true;
}
bool MonteCarloIpoptProblem::get_starting_point(
    Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z,
    Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda,
    Ipopt::Number* lambda) {
  if (not init_x) {
    throw std::runtime_error("init_x must be true");
  }
  if (init_z) {
    throw std::runtime_error("init_z must be false");
  }
  if (init_lambda) {
    throw std::runtime_error("init_lambda must be false");
  }

  // TODO: Use old initialization
  x[0] = (trajectory_param_min_ + trajectory_param_max_) / 2.0;
  // fmt::print("Initializing x to {}\n", x[0]);
  return true;
}

bool MonteCarloIpoptProblem::eval_f(Ipopt::Index n, const Ipopt::Number* x,
                                    bool new_x, Ipopt::Number& obj_value) {
  Recompute(x[0]);
  obj_value = last_cost_values_.cost_;
  return true;
}

bool MonteCarloIpoptProblem::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x,
                                         bool new_x, Ipopt::Number* grad_f) {
  Recompute(x[0]);
  grad_f[0] = last_cost_values_.cost_gradient_wrt_param_;
  return true;
}

bool MonteCarloIpoptProblem::eval_g(Ipopt::Index n, const Ipopt::Number* x,
                                    bool new_x, Ipopt::Index m,
                                    Ipopt::Number* g) {
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
    g[num_classical_constraints_] = last_risk_constraint_values_;
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

bool MonteCarloIpoptProblem::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x,
                                        bool new_x, Ipopt::Index m,
                                        Ipopt::Index nele_jac,
                                        Ipopt::Index* iRow, Ipopt::Index* jCol,
                                        Ipopt::Number* values) {
  if (values == NULL) {
    // Return structure of the jacobian of the constraints
    // element at i,j: grad_{x_j} g_{i}(x)
    for (Ipopt::Index r = 0; r < m; ++r) {
      for (Ipopt::Index c = 0; c < n; ++c) {
        Ipopt::Index idx = (r * n) + c;
        iRow[idx] = r;
        jCol[idx] = c;
      }
    }
  } else {
    throw std::runtime_error("Not providing jacobian ");
    //      Recompute(x[0]);
    //      if (classical_constraint_enabled_) {
    //        classical_constraints_.GetGradient(x[0], values,
    //                                           num_classical_constraints_);
    //      }
    //      if (risk_constraint_enabled_) {
    //        values[num_classical_constraints_] =
    //        last_risk_constraint_values_.d_risk_;
    //      }
  }
  return true;
}

bool MonteCarloIpoptProblem::eval_h(Ipopt::Index n, const Ipopt::Number* x,
                                    bool new_x, Ipopt::Number obj_factor,
                                    Ipopt::Index m, const Ipopt::Number* lambda,
                                    bool new_lambda, Ipopt::Index nele_hess,
                                    Ipopt::Index* iRow, Ipopt::Index* jCol,
                                    Ipopt::Number* values) {
  // return Ipopt::TNLP::eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda,
  // nele_hess, iRow, jCol, values);
  // if ((not risk_constraint_enabled_) or (classical_constraint_enabled_)) {
  //   return false;
  // }
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // element at 1,1: grad^2_{x1,x1} L(x,lambda)
    // TODO verify
    iRow[0] = 0;
    jCol[0] = 0;
  } else {
    throw std::runtime_error("Not providing hessian");
    // TODO classical constraints?
    // return the values
    // element at 1,1: grad^2_{x1,x1} L(x,lambda)
    // TODO verify

    //    Recompute(x[0]);
    //    values[0] = obj_factor * last_cost_values_.cost_hessian_wrt_param_;
    //    // TODO is this correct?
    //    if (risk_constraint_enabled_) {
    //      // if (risk_constraint_enabled_ and (not
    //      classical_constraint_enabled_)) { values[0] += lambda[0] *
    //      last_risk_constraint_values_.d2_risk_;
    //    }
  }
  return true;
}

void MonteCarloIpoptProblem::finalize_solution(
    Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
    const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
    const Ipopt::Number* g, const Ipopt::Number* lambda,
    Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data,
    Ipopt::IpoptCalculatedQuantities* ip_cq) {
  std::string out_str =
      "Finalized at " + std::to_string(x[0]) +
      " [risk time: " + std::to_string(timings_.risk_evals_seconds_) + " in " +
      std::to_string(timings_.num_risk_evals_) + " iters]\n";
  fmt::print("[IPOPT TNLP] Finalizing solution {} [status: {}]\n", x[0],
             ::roahm::ToString(status));
  final_return_status_ = status;
  // TODO or status == Ipopt::SolverReturn::FEASIBLE_POINT_FOUND?
  if (status == Ipopt::SolverReturn::SUCCESS) {
    fmt::print("STATUS WAS SUCCESS");
    found_feasible_internal_ = true;
    final_param_ = x[0];
    final_cost_ = obj_value;
    final_location_ = waypoint_cost_.GetPointAt(x[0]);
  }
}

bool MonteCarloIpoptProblem::FoundFeasible() const {
  return found_feasible_internal_;
}

std::optional<Ipopt::SolverReturn>
MonteCarloIpoptProblem::GetSolverReturnStatus() const {
  return final_return_status_;
}

std::optional<double> MonteCarloIpoptProblem::GetFeasibleCost() const {
  return final_cost_;
}

std::optional<double> MonteCarloIpoptProblem::GetFeasibleParam() const {
  return final_param_;
}
std::optional<PointXYH> MonteCarloIpoptProblem::GetFeasibleLocation() const {
  return final_location_;
}

WaypointLocalMirrorTakenIntoAccount
MonteCarloIpoptProblem::GetWaypointHeuristicAdjusted() const {
  return waypoint_heuristic_adjusted_;
}

WaypointLocalMirrorTakenIntoAccount
MonteCarloIpoptProblem::GetWaypointLocalMirrorAccounted() const {
  return waypoint_local_mirror_accounted_;
}

} // namespace roahm::monte_carlo_ipopt
