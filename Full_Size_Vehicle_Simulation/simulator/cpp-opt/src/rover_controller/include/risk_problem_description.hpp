#ifndef ROAHM_RISK_PROBLEM_DESCRIPTION_HPP_
#define ROAHM_RISK_PROBLEM_DESCRIPTION_HPP_

#include "fl_zono_obs_set.hpp"
#include "frs_select_info.hpp"
#include "frs_total.hpp"
#include "mu_sigma_multi.hpp"
#include "point_xyh.hpp"
#include "rover_state.hpp"
#include "unit_conversion.hpp"
#include "waypoint.hpp"
#include <optional>
#include <vector>

namespace roahm {

struct MaybeRisky {
  ::roahm::DynObs dyn_obs_;
  MuSigmaMulti mu_sigmas_;

  MaybeRisky(const ::roahm::DynObs& dyn_obs, const MuSigmaMulti& mu_sigmas)
      : dyn_obs_{dyn_obs}, mu_sigmas_{mu_sigmas} {}
};

struct RiskProblemDescription {
  /// @brief The initial rover state
  ::roahm::RoverState rover_state_;

  /// @brief Obstacles that will either use the REFINE or RISK constraints,
  /// depending on relative risk
  std::vector<MaybeRisky> maybe_risky_obs_;

  /// @brief Obstacles that we will ONLY use the RISK constraints on (e.g. on
  /// hardware)
  std::vector<MuSigmaMulti> always_risky_obs_;

  /// @brief Obstacles that we will ALWAYS use the REFINE constraints on (e.g.
  /// world bounds)
  ::roahm::FlZonoObsSet always_non_risky_obs_;

  /// @brief The desired waypoint for cost function
  ::roahm::WaypointGlobalNoMirror waypoint_global_no_mirror_;
};

struct RiskRtdPlanningOutputs {
  bool ran_;
  bool solve_succeeded_;
  bool found_feasible_;
  FrsSelectInfo sel_inf_{ManuType::kNone, -1, -1, -1, false};
  std::optional<double> final_param_;
  std::optional<double> final_cost_;
  std::optional<double> u0_;
  std::optional<PointXYH> final_location_;
  WaypointLocalMirrorTakenIntoAccount waypoint_with_heuristic_;
  WaypointLocalMirrorTakenIntoAccount waypoint_local_frame_;
};

struct OptimizationInputs {
  /// @brief Obstacle PDFs to generate risk constraints on.
  /// These should be in local frame, with mirroring taken into account.
  std::vector<std::vector<double>> mu_sigmas_;

  /// @brief Obstacles that have deterministic constraints.
  /// These should be in local frame, with mirroring taken into account.
  ::roahm::FlZonoObsSet fl_zono_obs_;

  /// @brief The desired waypoint in the local frame, with
  /// mirroring taken into account.
  ::roahm::WaypointLocalMirrorTakenIntoAccount waypoint_local_mirror_accounted_;

  ::roahm::WaypointLocalMirrorTakenIntoAccount waypoint_heuristic_adjusted_;
  /// @brief The starting speed of the ego vehicle, in meters per second.
  double u0_meters_per_second_;

  bool use_fixed_risk_threshold_;
  double risk_threshold_;
  bool use_left_cost_function_;

  /// @brief The number of risk constraints that will be present in the
  /// optimization problem.
  /// @return The number of risk constraints that will be present in the
  /// optimization problem.
  [[nodiscard]] int NumRiskConstraints() const noexcept {
    // TODO verify
    return (mu_sigmas_.empty()) ? 0 : 1;
  }

  /// @brief Returns true if there are any risk constraints.
  /// @return true if there are any risk constraints.
  [[nodiscard]] bool HasRiskConstraints() const noexcept {
    return NumRiskConstraints() > 0;
  }

  /// @brief Returns true if there are any classical/deterministic constraints.
  /// @return true if there are any classical/deterministic constraints.
  [[nodiscard]] bool HasClassicalConstraints() const noexcept {
    return fl_zono_obs_.GetNumObs() > 0;
  }

  /// @brief Returns the starting speed of the ego vehicle in meters per second.
  /// @return the starting speed of the ego vehicle in meters per second.
  [[nodiscard]] double GetU0MetersPerSecond() const noexcept {
    return u0_meters_per_second_;
  }
};

inline static ::roahm::OptimizationInputs
GetOptInputs(const ::roahm::FrsTotal& frs_total,
             const ::roahm::FrsSelectInfo& select_info,
             const ::roahm::RiskProblemDescription& base_inputs,
	     const bool use_waypoint_modification_heuristic,
	     const bool use_fixed_risk_threshold,
	     const double risk_threshold,
	     const bool use_left_cost_function) {
  const auto init_state = base_inputs.rover_state_;
  const double u0 = init_state.u_;
  const bool mirror = select_info.mirror_;
  const auto& vehrs = frs_total.GetVehrs(select_info);

  ::roahm::OptimizationInputs ret;
  ret.use_fixed_risk_threshold_ = use_fixed_risk_threshold;
  ret.risk_threshold_ = risk_threshold;
  ret.use_left_cost_function_ = use_left_cost_function;
  ret.u0_meters_per_second_ = u0;

  ret.waypoint_local_mirror_accounted_ =
      base_inputs.waypoint_global_no_mirror_
          .RelativeToEgoFrameNoMirror(base_inputs.rover_state_.GetXYH())
          .TakeMirrorIntoAccount(mirror);

  ret.waypoint_heuristic_adjusted_ = use_waypoint_modification_heuristic ? 
      ret.waypoint_local_mirror_accounted_.HeuristicAdjust(vehrs.GetManuType(),
                                                           u0) : ret.waypoint_local_mirror_accounted_;

  for (const auto& always_risky_mu_sigma_multi :
       base_inputs.always_risky_obs_) {
    ret.mu_sigmas_.push_back(vehrs.ExtractMuSigma(always_risky_mu_sigma_multi,
                                                  init_state.GetXYH(), mirror));
  }

  // ret.fl_zono_obs_ = base_inputs.always_non_risky_obs_;
  ret.fl_zono_obs_ = {};
  for (int i = 0; i < base_inputs.always_non_risky_obs_.GetNumObs(); ++i) {
    ret.fl_zono_obs_.PushObs(
        base_inputs.always_non_risky_obs_.GetSingleDynObs(i).RelativeTo(
            init_state.GetXYH(), mirror));
  }
  // TODO REMOVE for (const auto& always_non_risky_dyn_obs :
  // TODO REMOVE      base_inputs.always_non_risky_obs_) {
  // TODO REMOVE   ret.fl_zono_obs_.PushObs(
  // TODO REMOVE       always_non_risky_dyn_obs.RelativeTo(init_state.GetXYH(),
  // mirror));
  // TODO REMOVE }

  {
    int maybe_risky_idx{0};
    for (const auto& maybe_risky : base_inputs.maybe_risky_obs_) {
      // TODO fix this
      const double obs_vel{maybe_risky.dyn_obs_.Velocity()};
      const double v_rel = std::abs(u0 - obs_vel);
      // CHALLEN_NOTE set this to +1000.0 for RISK, -1000.0 for REFINE
      constexpr double kRiskRelVelocityMaxMetersPerSecond =
          MilesPerHourToMetersPerSecond(1000.0);
      const bool is_risky = v_rel <= kRiskRelVelocityMaxMetersPerSecond;
      if (is_risky) {
        ret.mu_sigmas_.push_back(vehrs.ExtractMuSigma(
            maybe_risky.mu_sigmas_, init_state.GetXYH(), mirror));
      } else {
        const auto obs_rel =
            maybe_risky.dyn_obs_.RelativeTo(init_state.GetXYH(), mirror);
        // fmt::print("Non Risky: [x: {}] [y: {}] [h: {}] [vel: {}] [len: {}] "
        //            "[width: {}]\n",
        //            obs_rel.CenterX0(), obs_rel.CenterY0(),
        //            obs_rel.HeadingRad(), obs_rel.Velocity(),
        //            obs_rel.Length(), obs_rel.Width());
        ret.fl_zono_obs_.PushObs(obs_rel);
      }
      // fmt::print("ChRisky [{}, v_rel: {}, obs vel: {}]: {}\n",
      // maybe_risky_idx,
      //            v_rel, obs_vel, is_risky);
      ++maybe_risky_idx;
    }
  }
  return ret;
}

} // namespace roahm
#endif // ROAHM_RISK_PROBLEM_DESCRIPTION_HPP_
