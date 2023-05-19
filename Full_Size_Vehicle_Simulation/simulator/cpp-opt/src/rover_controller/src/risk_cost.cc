#include "risk_cost.hpp"

namespace roahm::risk_cost {

inline double Square(double val) { return val * val; }

RiskCost::RiskCost(const double traj_param_min, const double traj_param_max)
    : traj_param_min_{traj_param_min}, traj_param_max_{traj_param_max} {}

RiskCostValues RiskCost::EvaluateAt(double traj_param) const {
  const double traj_param_range = traj_param_max_ - traj_param_min_;
  const double b = (0.75 * traj_param_range) + traj_param_min_;
  RiskCostValues values;
  values.cost_ = Square(traj_param - b);
  values.cost_gradient_ = 2.0 * (traj_param - b);
  values.cost_hessian_ = 2.0;
  return values;
}
} // namespace roahm::risk_cost