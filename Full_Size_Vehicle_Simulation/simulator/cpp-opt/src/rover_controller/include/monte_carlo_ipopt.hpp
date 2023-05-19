#ifndef ROAHM_MONTE_CARLO_IPOPT_HPP_
#define ROAHM_MONTE_CARLO_IPOPT_HPP_
#include "IpAlgTypes.hpp" // for SolverReturn
#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"  // for IpoptCalculatedQuantities, IpoptData, TNLP
#include "IpTypes.hpp" // for Index, Number, Ipopt
#include "comparison_methods.hpp"
#include "risk_problem_description.hpp"
#include "rover_state.hpp"
namespace roahm::monte_carlo_ipopt {

class MonteCarloIpoptProblem : public Ipopt::TNLP {
private:
  /// Convenience to avoid fully qualifying type
  using Number = Ipopt::Number;
  /// Convenience to avoid fully qualifying type
  using Index = Ipopt::Index;

  // Risk constraints
  bool risk_constraint_enabled_;
  int num_risk_constraints_;
  // Decision variable bounds
  double trajectory_param_min_;
  double trajectory_param_max_;
  // Classical FL Zono constraints
  ::roahm::fl_zono_constraint::FlZonoConstraint classical_constraints_;
  int num_classical_constraints_;
  bool classical_constraint_enabled_;
  // Cached values
  double last_trajectory_param_;
  ::roahm::WaypointCostRet last_cost_values_;
  // Cost function
  ::roahm::WaypointCost waypoint_cost_;
  struct OptTimings {
    int num_risk_evals_;
    double risk_evals_seconds_;
  };
  OptTimings timings_;
  std::mt19937 gen_;

  std::optional<Ipopt::SolverReturn> final_return_status_;
  std::optional<double> final_param_;
  std::optional<double> final_cost_;
  std::optional<PointXYH> final_location_;
  ::roahm::WaypointLocalMirrorTakenIntoAccount waypoint_local_mirror_accounted_;
  ::roahm::WaypointLocalMirrorTakenIntoAccount waypoint_heuristic_adjusted_;
  bool found_feasible_internal_;
  const Vehrs unsliced_vehrs_;
  const RoverState z0_;
  const ::roahm::risk_comparisons::EnvFootprints footprints_;
  const std::vector<std::vector<double>> all_mu_sigmas_;
  const std::int64_t num_trajectory_samples_;
  double last_risk_constraint_values_;

  void Recompute(const double trajectory_param);

public:
  MonteCarloIpoptProblem(const Vehrs& vehrs,
                         const std::array<double, 3>& u0v0r0_slice_beta,
                         const OptimizationInputs& opt_inputs,
                         const std::size_t seed_val,
                         const Vehrs& unsliced_vehrs, const RoverState& z0,
                         const ::roahm::risk_comparisons::EnvFootprints footprints,
                         const int num_trajectory_samples);

  /// Sets information about the problem
  /// \param n out param to store dimensionality of the decision variable
  /// \param m out param to store the number of constraints
  /// \param nnz_jac_g out param to store number of nonzero Jacobian entries
  /// \param nnz_h_lag out param to store number of nonzero Hessian entries
  /// \param index_style whether to 0 (C style) or 1-indexed (Fortran style)
  /// \return true iff successful
  bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag,
                    IndexStyleEnum& index_style) override;

  /// Sets bounds on variables and constraints
  /// \param n dimensionality of the decision variable, x
  /// \param x_l lower bounds on each dimension of the decision variable
  /// \param x_u upper bounds on each dimension of the decision variable
  /// \param m number of constraints
  /// \param g_l lower bounds for each of the constraints
  /// \param g_u upper bounds for each of the constraints
  /// \return true iff successful
  bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l,
                       Number* g_u) override;

  /// Gets the starting point for the optimization
  /// \param n dimensionality of the decision variable, x
  /// \param init_x if true, an initial guess must be provided for \p x
  /// \param x out param for providing initial guess for the decision variable
  /// \param init_z if true, initial guess for bound multipliers \p z_l and
  /// \p z_u must be provided
  /// \param z_L out param for initial value for lower bound for bound
  /// multipliers \param z_U out param for initial value for upper bound for
  /// bound multipliers \param m number of constraints \param init_lambda if
  /// true, must provide initial value for the constraint multipliers \p lambda
  /// \param lambda out param for initial value for constraint multipliers
  /// \return true iff successful
  bool get_starting_point(Index n, bool init_x, Number* x, bool init_z,
                          Number* z_L, Number* z_U, Index m, bool init_lambda,
                          Number* lambda) override;

  /// Evaluates the cost function to compute the objective value
  /// \param n the dimensionality of the decision variable, \p x
  /// \param x the evaluation point to call \f$ f(x) \f$
  /// \param new_x false if any eval_* function has been called with the same
  /// value of \p x, true otherwise
  /// \param obj_value out param to store the value of the objective function
  /// \return true iff successful
  bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) override;

  /// Computes the gradient of the objective function
  /// \param n dimensionality of the decision variable, \p x
  /// \param x the value of the decision variable to evaluate the gradient at
  /// \param new_x false if any eval_* method was called with the same value
  /// in \p x, true otherwise
  /// \param grad_f out param to store the values of the gradient of the
  /// objective function
  /// \return true iff successful
  bool eval_grad_f(Index n, const Number* x, bool new_x,
                   Number* grad_f) override;

  /// Computes the constraint values
  /// \param n the dimensionality of the decision variable, \p x
  /// \param x the decision variable
  /// \param new_x false if any eval_* was called with the same value in \p x
  /// \param m the number of constraints
  /// \param g out param to store the constraint values
  /// \return true iff successful
  bool eval_g(Index n, const Number* x, bool new_x, Index m,
              Number* g) override;

  /// Evaluates the jacobian of the gradient, and sets the structure (location
  /// of nonzero entries) if it is the first evaluation
  /// \param n the dimensionality of the decision variable, \p x
  /// \param x the decision variable
  /// \param new_x false if any eval_* was called with the same value in \p x
  /// \param m the number of constraints
  /// \param nele_jac the number of nonzero elements in the Jacobian
  /// \param iRow out param. First call: array of length \p nele_jac to store
  /// the row indices of entries in the Jacobian of the constraints.
  /// Later calls: NULL
  /// \param jCol out param. First call: array of length \p nele_jac to store
  /// the column indices of entries in the Jacobian of the constraints.
  /// Later calls: NULL
  /// \param values out param. First call: NULL. Later calls: array of length
  /// \p nele_jac to store the values of the entries in the Jacobian of the
  /// constraints
  /// \return true iff successful
  bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac,
                  Index* iRow, Index* jCol, Number* values) override;

  /// First call: returns the structure of the hessian of the lagrangian.
  /// Secondary calls: evaluates the hessian of the lagrangian.
  /// \param n the dimensionality of the decision variable, \p x
  /// \param x the decision variable
  /// \param new_x false if any eval_* has been called with the same \p x
  /// \param obj_factor the factor in front of the objective term in the
  /// Hessian
  /// \param m the number of constraints
  /// \param lambda the constraint multipliers to evaluate the Hessian at
  /// \param new_lambda false if any eval_* method was called with the same
  /// values in \p lambda
  /// \param nele_hess the number of nonzero elements in the Hessian
  /// \param iRow out param. First call: array of length \p nele_hess to store
  /// row indices of entries in the Hessian. Later calls: NULL
  /// \param jCol out param. First call: array of length \p nele_hess to store
  /// column indices of entries in the Hessian. Later calls: NULL
  /// \param values out param. First call: NULL. Later calls: array of length
  /// \p nele_hess to store values of entries in the Hessian
  /// \return true iff successful
  bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m,
              const Number* lambda, bool new_lambda, Index nele_hess,
              Index* iRow, Index* jCol, Number* values) override;

  /// Called when IPOPT if finished to handle the solution or failure cases
  /// \param status provides the status of the algorithm
  /// \param n the dimensionality of the decision variable, \p x
  /// \param x the final decision variable
  /// \param z_L the final lower bound multipliers
  /// \param z_U the final upper bound multipliers
  /// \param m the number of constraints
  /// \param g the final constraint function values
  /// \param lambda the final constraint multipliers
  /// \param obj_value the final objective function value
  /// \param ip_data extra data
  /// \param ip_cq extra data
  void finalize_solution(Ipopt::SolverReturn status, Index n, const Number* x,
                         const Number* z_L, const Number* z_U, Index m,
                         const Number* g, const Number* lambda,
                         Number obj_value, const Ipopt::IpoptData* ip_data,
                         Ipopt::IpoptCalculatedQuantities* ip_cq) override;

  /// Returns true iff a feasible solution has been found at any point
  /// \return true iff a feasible solution has been found at any point
  bool FoundFeasible() const;

  std::optional<Ipopt::SolverReturn> GetSolverReturnStatus() const;

  /// Returns the cost associated with the minimum-cost feasible solution, if
  /// a feasible solution has been found at any point.
  /// \return if a feasible solution has been found at any point, returns the
  /// cost associated with the minimum-cost feasible solution, otherwise there
  /// are no guarantees on return value.
  std::optional<double> GetFeasibleCost() const;

  /// Returns the solution associated with the minimum-cost feasible solution,
  /// if a feasible solution has been found at any point.
  /// \return if a feasible solution has been found at any point, returns the
  /// solution associated with the minimum-cost feasible solution that has been
  /// found, otherwise there are no guarantees on the return value.
  std::optional<double> GetFeasibleParam() const;

  std::optional<PointXYH> GetFeasibleLocation() const;

  Waypoint<true, true> GetWaypointHeuristicAdjusted() const;
  Waypoint<true, true> GetWaypointLocalMirrorAccounted() const;

  ~MonteCarloIpoptProblem() override = default;
  MonteCarloIpoptProblem() = delete;
  MonteCarloIpoptProblem(const MonteCarloIpoptProblem&) = delete;
  MonteCarloIpoptProblem& operator=(const MonteCarloIpoptProblem&) = delete;
};
} // namespace roahm::monte_carlo_ipopt
#endif // ROAHM_MONTE_CARLO_IPOPT_HPP_