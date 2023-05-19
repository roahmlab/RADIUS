#ifndef ROAHM_WAYPOINT_COST_HPP_
#define ROAHM_WAYPOINT_COST_HPP_

#include "vehrs.hpp"
#include "waypoint.hpp"

namespace roahm {

/// @brief Contains the outputs from waypoint cost evaluation: the cost,
/// cost gradient, and cost hessian, with respect to the parameter input
struct WaypointCostRet {
  /// @brief The cost at the evaluated trajectory parameter.
  double cost_;

  /// @brief The cost gradient at the evaluation point, with respect to the
  /// trajectory parameter.
  double cost_gradient_wrt_param_;

  /// @brief The cost hessian at the evaluation point, with respect to the
  /// trajectory parameter.
  double cost_hessian_wrt_param_;
};

struct WaypointCost {
private:
  [[nodiscard]] double GetLambda(const double param) const;

public:
  /// @brief The goal waypoint.
  const WaypointLocalMirrorTakenIntoAccount des_wp_;

  /// @brief The center value along the sliceable trajectory parameter's
  /// dimension.
  const double sliceable_center_;

  /// @brief The generator value along the sliceable trajectory parameter's
  /// dimension.
  const double sliceable_generator_val_;

  /// @brief The x center value of the target zonotope.
  const double x_center_;

  /// @brief The y center value of the target zonotope.
  const double y_center_;

  /// @brief The heading center value of the target zonotope.
  const double h_center_;

  /// @brief The x generator value of the target zonotope on the trajectory
  /// parameter's sliceable generator.
  const double sliceable_x_generator_;

  /// @brief The y generator value of the target zonotope on the trajectory
  /// parameter's sliceable generator.
  const double sliceable_y_generator_;

  /// @brief The h generator value of the target zonotope on the trajectory
  /// parameter's sliceable generator.
  const double sliceable_h_generator_;

  /// @brief The weighting associated with the x portion of the cost function
  const double x_weight_;

  /// @brief The weighting associated with the y portion of the cost function
  const double y_weight_;

  /// @brief The weighting associated with the heading portion of the cost
  /// function
  const double h_weight_;

  /// @brief The epsilon to add to the x portion's square root, to avoid
  /// division by zero. Larger values are more inaccurate, but very small values
  /// may cause numerical issues and explode the gradient/hessian.
  const double x_sqrt_epsilon_;

  /// @brief The epsilon to add to the y portion's square root, to avoid
  /// division by zero. Larger values are more inaccurate, but very small values
  /// may cause numerical issues and explode the gradient/hessian.
  const double y_sqrt_epsilon_;

  /// @brief The epsilon to add to the heading portion's square root, to avoid
  /// division by zero. Larger values are more inaccurate, but very small values
  /// may cause numerical issues and explode the gradient/hessian.
  const double h_sqrt_epsilon_;

  /// @brief Constructor.

  WaypointCost(const WaypointLocalMirrorTakenIntoAccount& des_wp,
               const double sliceable_center,
               const double sliceable_generator_val, const double x_center,
               const double y_center, const double h_center,
               const double sliceable_x_generator,
               const double sliceable_y_generator,
               const double sliceable_h_generator, const double x_weight = 3.0,
               const double y_weight = 10.0, const double h_weight = 10.0,
               const double x_sqrt_epsilon = 1.0e-8,
               const double y_sqrt_epsilon = 1.0e-8,
               const double h_sqrt_epsilon = 1.0e-8);

  /// @brief Evaluates the cost at a given trajectory parameter
  /// @param param the trajectory parameter to evaluate the cost at. This value
  /// should be within the slicing range.
  /// @return the cost, gradient, and hessian at the evaluation point.
  WaypointCostRet EvaluateAt(const double param) const;

  PointXYH GetPointAt(const double param) const;
};

// TODO move to member fcn?
inline static double GetLookupTimeSeconds(const ManuType manu_type) {
  if (IsLan(manu_type)) {
    return 6.0;
  }
  return 3.0;
}

inline static WaypointCost WaypointCostFromVehrsAndWp(
    const ::roahm::Vehrs& vehrs,
    const WaypointLocalMirrorTakenIntoAccount& des_wp_modified) {
  const auto time_ahead = GetLookupTimeSeconds(vehrs.GetManuType());
  const auto zono_slice_xyh_cen_gen =
      vehrs.GetSliceableCenterAndGensNearstToTime(time_ahead);
  const auto zono_center_gen_traj_param = vehrs.GetCenterGenTrajParam();
  return ::roahm::WaypointCost{des_wp_modified,
                               zono_center_gen_traj_param.first,
                               zono_center_gen_traj_param.second,
                               zono_slice_xyh_cen_gen.first.x_,
                               zono_slice_xyh_cen_gen.first.y_,
                               zono_slice_xyh_cen_gen.first.h_,
                               zono_slice_xyh_cen_gen.second.x_,
                               zono_slice_xyh_cen_gen.second.y_,
                               zono_slice_xyh_cen_gen.second.h_};
}

} // namespace roahm
#endif // ROAHM_WAYPOINT_COST_HPP_