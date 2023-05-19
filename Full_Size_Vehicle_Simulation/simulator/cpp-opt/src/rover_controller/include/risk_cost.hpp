#ifndef ROAHM_RISK_COST_HPP_
#define ROAHM_RISK_COST_HPP_
namespace roahm::risk_cost {
struct RiskCostValues {
  double cost_;
  double cost_gradient_;
  double cost_hessian_;
};

struct RiskCost {
  const double traj_param_min_;
  const double traj_param_max_;
  RiskCost(const double traj_param_min, const double traj_param_max);

  RiskCostValues EvaluateAt(double traj_param) const;
};
} // namespace roahm::risk_cost
#endif