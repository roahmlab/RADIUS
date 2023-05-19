#include "waypoint_cost.hpp"

#include "waypoint.hpp"

namespace roahm {
WaypointCost::WaypointCost(
    const WaypointLocalMirrorTakenIntoAccount& des_wp,
    const double sliceable_center, const double sliceable_generator_val,
    const double x_center, const double y_center, const double h_center,
    const double sliceable_x_generator, const double sliceable_y_generator,
    const double sliceable_h_generator, const double x_weight,
    const double y_weight, const double h_weight, const double x_sqrt_epsilon,
    const double y_sqrt_epsilon, const double h_sqrt_epsilon)
    : des_wp_{des_wp}, sliceable_center_{sliceable_center},
      sliceable_generator_val_{sliceable_generator_val}, x_center_{x_center},
      y_center_{y_center}, h_center_{h_center},
      sliceable_x_generator_{sliceable_x_generator},
      sliceable_y_generator_{sliceable_y_generator},
      sliceable_h_generator_{sliceable_h_generator}, x_weight_{x_weight},
      y_weight_{y_weight}, h_weight_{h_weight}, x_sqrt_epsilon_{x_sqrt_epsilon},
      y_sqrt_epsilon_{y_sqrt_epsilon}, h_sqrt_epsilon_{h_sqrt_epsilon} {}

[[nodiscard]] double WaypointCost::GetLambda(const double param) const {
  return (param - sliceable_center_) / sliceable_generator_val_;
}
PointXYH WaypointCost::GetPointAt(const double param) const {
  const double lambda{GetLambda(param)};
  const double x_val{x_center_ + (lambda * sliceable_x_generator_)};
  const double y_val{y_center_ + (lambda * sliceable_y_generator_)};
  const double h_val{h_center_ + (lambda * sliceable_h_generator_)};
  return PointXYH{x_val, y_val, h_val};
}
WaypointCostRet WaypointCost::EvaluateAt(const double param) const {
  // Cost as a function of lamdba:
  //
  // C(lambda)
  // = C_xy(lambda) + C_h(lambda)
  // = W_xy * sqrt((A-B*lambda)^2 + (C - D*lambda)^2 + K)
  //  + W_h sqrt((E - F*lambda)^2 + K2)

  const double lambda{GetLambda(param)};
  const double A = des_wp_.x_ - x_center_;
  const double B = sliceable_x_generator_;
  const double C = des_wp_.y_ - y_center_;
  const double D = sliceable_y_generator_;
  const double E = des_wp_.h_ - h_center_;
  const double F = sliceable_h_generator_;
  const double K0 = x_sqrt_epsilon_;
  const double K1 = x_sqrt_epsilon_;
  const double K2 = h_sqrt_epsilon_;

  // Compute error on each state
  const double err_x = A - lambda * B;
  const double err_y = C - lambda * D;
  const double err_h = E - lambda * F;

  // Square of errors
  const double err_x2 = err_x * err_x;
  const double err_y2 = err_y * err_y;
  const double err_h2 = err_h * err_h;

  // const double unscaled_cost_xy_unsq = err_x2 + err_y2 + K;
  // const double unscaled_cost_xy = std::sqrt(unscaled_cost_xy_unsq);
  const double unscaled_cost_x = std::sqrt(err_x2 + K0);
  const double unscaled_cost_y = std::sqrt(err_y2 + K1);
  const double unscaled_cost_h = std::sqrt(err_h2 + K2);
  const double cost_eval = (x_weight_ * unscaled_cost_x) +
                           (y_weight_ * unscaled_cost_y) +
                           (h_weight_ * unscaled_cost_h);

  // Derivative of lambda with respect to the sliceable parameter
  const double dl_dk = 1.0 / sliceable_generator_val_;

  // Derviatives of the x, y, h cost with respect to lambda
  const double d_unscaled_cost_x_dl =
      (B * (B * lambda - A)) / (std::sqrt(err_x2 + K0));

  const double d_unscaled_cost_y_dl =
      (D * (D * lambda - C)) / (std::sqrt(err_y2 + K1));

  const double d_unscaled_cost_h_dl =
      (F * (F * lambda - E)) / (std::sqrt(err_h2 + K2));

  // Derivative of the cost with respect to lambda
  const double d_cost_dl = x_weight_ * d_unscaled_cost_x_dl +
                           y_weight_ * d_unscaled_cost_y_dl +
                           h_weight_ * d_unscaled_cost_h_dl;

  // Derivative of the cost with respect to the sliceable parameter
  const double d_cost_dk = d_cost_dl * dl_dk;

  const double B2 = B * B;
  const double D2 = D * D;
  const double F2 = F * F;

  const double d_cost_x2_dl2 =
      (B2 * K0) / (unscaled_cost_x * unscaled_cost_x * unscaled_cost_x);
  const double d_cost_y2_dl2 =
      (D2 * K1) / (unscaled_cost_y * unscaled_cost_y * unscaled_cost_y);
  const double d_cost_h2_dl2 =
      (F2 * K2) / (unscaled_cost_h * unscaled_cost_h * unscaled_cost_h);
  const double d_cost2_dl2 = x_weight_ * d_cost_x2_dl2 +
                             y_weight_ * d_cost_y2_dl2 +
                             h_weight_ * d_cost_h2_dl2;
  const double d_cost2_dk2 = d_cost2_dl2 * dl_dk * dl_dk;

  return {cost_eval, d_cost_dk, d_cost2_dk2};
}
} // namespace roahm