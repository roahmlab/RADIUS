#ifndef ROAHM_MONTE_CARLO_OPT_E_HPP_
#define ROAHM_MONTE_CARLO_OPT_E_HPP_
#include "fl_zono_constraint.hpp"
#include "fl_zono_obs_set.hpp"
#include "frs_loader.hpp"
#include "frs_mega.hpp"
#include "frs_select_info.hpp"
#include "frs_total.hpp"
#include "individual_mu_sigma.hpp"
#include "mu_sigma_multi.hpp"
#include "point_xyh.hpp"
#include "risk_constraint.hpp"
#include "risk_problem_description.hpp"
#include "rover_state.hpp"
#include "timing_util.hpp"
#include "waypoint.hpp"
#include "waypoint_cost.hpp"
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <fmt/format.h>
#include <random>
#include <stdexcept>

namespace roahm::risk_comparisons {

enum class CollisionStatus { kCollisionOccurred, kNoCollisionOccurred };

struct Footprint2d {
  double length_;
  double width_;
};

struct EnvFootprints {
  Footprint2d ego_footprint_;
  Footprint2d obs_footprint_;

  static inline EnvFootprints Default() {
    return {.ego_footprint_{2.4, 1.1}, .obs_footprint_{2.4, 1.1}};
  }
};

class PolytopeHalfspaceConstraintMats {
public:
  Eigen::Matrix<double, Eigen::Dynamic, 2> a_mat_;
  Eigen::VectorXd b_mat_;
};

inline static PolytopeHalfspaceConstraintMats
PolytopeHalfspace(Eigen::Vector2d center,
                  Eigen::Matrix<double, 2, Eigen::Dynamic> gens) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> c;
  c.conservativeResize(2, gens.cols());
  c.row(0) = -gens.row(1);
  c.row(1) = gens.row(0);
  c.colwise().normalize();
  c.transposeInPlace();
  Eigen::Matrix<double, Eigen::Dynamic, 1> delta_d =
      ((((c * gens).transpose()).array().abs()).colwise().sum()).transpose();
  Eigen::Matrix<double, Eigen::Dynamic, 1> d = c * center;

  PolytopeHalfspaceConstraintMats cons;
  auto& a_mat{cons.a_mat_};
  a_mat.conservativeResize(gens.cols() * 2, 2);
  a_mat << c, -c;
  auto& b_mat{cons.b_mat_};
  b_mat.conservativeResize(gens.cols() * 2, 1);
  auto b_mat1 = delta_d + d;
  auto b_mat2 = delta_d - d;
  b_mat << b_mat1, b_mat2;
  return cons;
}

inline static CollisionStatus CheckIntersection(Eigen::Vector2d ego_xy,
                                                Eigen::Vector2d obs_xy,
                                                Eigen::Matrix2d ego_gens,
                                                Eigen::Matrix2d obs_gens) {
  const Eigen::Vector2d delta_xy{obs_xy - ego_xy};
  Eigen::Matrix<double, 2, 4> gens;
  gens << ego_gens, obs_gens;
  const auto cons{PolytopeHalfspace(delta_xy, gens)};
  return (0.0 <= cons.b_mat_.minCoeff())
             ? CollisionStatus::kCollisionOccurred
             : CollisionStatus::kNoCollisionOccurred;
}

inline static CollisionStatus
CheckIntersectionArbitrary(Eigen::Vector2d ego_xy, Eigen::Vector2d obs_xy,
                           Eigen::Matrix2Xd ego_gens, Eigen::Matrix2Xd obs_gens,
                           std::optional<double> ego_buff_x,
                           std::optional<double> ego_buff_y) {
  const Eigen::Vector2d delta_xy{obs_xy - ego_xy};
  Eigen::Matrix2Xd gens;
  std::int64_t num_gens = ego_gens.cols() + obs_gens.cols() +
                          (ego_buff_x.has_value()) + (ego_buff_y.has_value());
  gens.conservativeResize(2, num_gens);
  std::int64_t last_col{0};
  gens.block(0, last_col, 2, ego_gens.cols()) = ego_gens;
  last_col += ego_gens.cols();
  if (ego_buff_x.has_value()) {
    Eigen::Vector2d ego_buff_x_gen{ego_buff_x.value(), 0.0};
    gens.col(last_col) = ego_buff_x_gen;
    ++last_col;
  }
  if (ego_buff_y.has_value()) {
    Eigen::Vector2d ego_buff_y_gen{0.0, ego_buff_y.value()};
    gens.col(last_col) = ego_buff_y_gen;
    ++last_col;
  }
  gens.block(0, last_col, 2, obs_gens.cols()) = obs_gens;
  const auto cons{PolytopeHalfspace(delta_xy, gens)};
  return (0.0 <= cons.b_mat_.minCoeff())
             ? CollisionStatus::kCollisionOccurred
             : CollisionStatus::kNoCollisionOccurred;
}

inline static Eigen::Matrix2d MakeEigMat2d(const double m00, const double m01,
                                           const double m10, const double m11) {
  Eigen::Matrix2d ret;
  ret << m00, m01, m10, m11;
  return ret;
}

inline static double SampleNormal(const double standard_deviation,
                                  std::mt19937& gen) {
  std::normal_distribution<double> d{0.0, standard_deviation};
  return d(gen);
}

inline static double SampleUniform(const double min_val, const double max_val,
                                   std::mt19937& gen) {
  std::uniform_real_distribution<double> d{min_val, max_val};
  return d(gen);
}

inline static Eigen::Vector2d
GetSample(const double mu_x, const double mu_y, const double sigma_xx,
          const double sigma_xy, const double sigma_yy, std::mt19937& gen) {

  const Eigen::Matrix2d covar_unrot{
      MakeEigMat2d(sigma_xx, sigma_xy, sigma_xy, sigma_yy)};
  Eigen::EigenSolver<Eigen::Matrix2d> eigen_solver;
  eigen_solver.compute(covar_unrot);
  const auto covar_eigen_vals{eigen_solver.eigenvalues().real()};
  const auto covar_eigen_vecs{eigen_solver.eigenvectors().real().eval()};
  if (covar_eigen_vals.x() < 0.0 or covar_eigen_vals.y() < 0.0) {
    throw std::runtime_error("Eigen values are negative");
  }
  const auto stdevs_rot{covar_eigen_vals.cwiseSqrt().eval()};
  const double delta_x_rot{SampleNormal(stdevs_rot.x(), gen)};
  const double delta_y_rot{SampleNormal(stdevs_rot.y(), gen)};
  const Eigen::Vector2d delta_xy_rot{delta_x_rot, delta_y_rot};
  const Eigen::Vector2d delta_xy_real =
      (covar_eigen_vecs * delta_xy_rot).eval();
  const Eigen::Vector2d mu_xy{mu_x, mu_y};
  const Eigen::Vector2d ret{mu_xy + delta_xy_real};
  return ret;
}

inline static Eigen::Vector2d SampleMuSigma(const IndividualMuSigma& mu_sigma,
                                            std::mt19937& gen) {
  return GetSample(mu_sigma.mu_x_, mu_sigma.mu_y_, mu_sigma.sigma_1_,
                   mu_sigma.sigma_2_, mu_sigma.sigma_4_, gen);
}

inline static std::vector<IndividualMuSigma>
MuSigmaIndividFromRawVec(const std::vector<double>& mu_sigma_vec) {
  if (mu_sigma_vec.size() % 5 != 0) {
    throw std::runtime_error("Raw mu sigma vec has size not divisible by 5");
  }
  std::vector<IndividualMuSigma> mu_sigma_individ_vec;
  for (int i = 0; i < mu_sigma_vec.size(); i += 5) {
    const double mu_x{mu_sigma_vec.at(i + 0)};
    const double mu_y{mu_sigma_vec.at(i + 1)};
    const double sigma_xx{mu_sigma_vec.at(i + 2)};
    const double sigma_xy{mu_sigma_vec.at(i + 3)};
    const double sigma_yy{mu_sigma_vec.at(i + 4)};
    mu_sigma_individ_vec.push_back({.mu_x_ = mu_x,
                                    .mu_y_ = mu_y,
                                    .sigma_1_ = sigma_xx,
                                    .sigma_2_ = sigma_xy,
                                    .sigma_4_ = sigma_yy});
  }
  return mu_sigma_individ_vec;
}

inline static Eigen::Matrix2Xd
SampleSingleTrajectory(const MuSigmaMulti& mu_sigmas,
                       const TimeMicroseconds traj_duration,
                       const TimeMicroseconds dt, std::mt19937& gen,
                       const PointXYH& local_frame, const bool mirror) {
  const auto mu_sig_global{
      MuSigmaIndividFromRawVec(mu_sigmas.GetMuSigmaVec(dt))};
  Eigen::Matrix2Xd ret;
  ret.conservativeResize(2, traj_duration / dt + 1);
  {
    std::int64_t i{0};
    for (TimeMicroseconds time_since_traj_start{
             TimeMicroseconds::FromMicroseconds(0)};
         time_since_traj_start <= traj_duration; time_since_traj_start += dt)
      ret.col(i++) = SampleMuSigma(
          mu_sig_global.at(i).RelativeTo(local_frame, mirror), gen);
  }
  return ret;
}

inline static Eigen::Matrix3Xd
GetEgoTrajectoryGlobalFrame(const Vehrs& sliced_vehrs,
                            const TimeMicroseconds dt, const PointXYH z0) {
  Eigen::Matrix3Xd ret;
  const double h0{z0.h_};
  const Eigen::Vector2d ego_xy_global{z0.x_, z0.y_};
  const Eigen::Matrix2d rot_h0{Eigen::Rotation2Dd(z0.h_).matrix()};
  const TimeMicroseconds traj_duration{sliced_vehrs.TrajectoryDuration()};
  ret.conservativeResize(3, traj_duration / dt + 1);
  {
    std::int64_t i{0};
    for (TimeMicroseconds duration_since_traj_start{
             TimeMicroseconds::FromMicroseconds(std::int64_t{0})};
         duration_since_traj_start <= traj_duration;
         duration_since_traj_start += dt) {
      const auto interp_center{
          sliced_vehrs.LerpCenter(duration_since_traj_start.ToSecondsDouble())};
      const Eigen::Vector2d zono_center_xy_wrt_frs{interp_center.x_,
                                                   interp_center.y_};
      const Eigen::Vector2d zono_center_xy_global{
          (rot_h0 * zono_center_xy_wrt_frs) + ego_xy_global};
      ret.col(i++) =
          Eigen::Vector3d{zono_center_xy_global.x(), zono_center_xy_global.y(),
                          h0 + interp_center.h_};
    }
  }
  return ret;
}

// return true IFF collision
inline static CollisionStatus
CollisionCheckTrajectory(const Eigen::Matrix3Xd& ego_trajectory_global,
                         const Eigen::Matrix2Xd& obs_trajectory_global,
                         const EnvFootprints& footprints) {
  if (ego_trajectory_global.cols() != obs_trajectory_global.cols()) {
    throw std::runtime_error(fmt::format(
        "Number of time instances in ego and obstacle trajectory not equal [{} "
        "!= {}]",
        ego_trajectory_global.cols(), obs_trajectory_global.cols()));
  }

  const Eigen::Matrix<double, 2, 2> obs_generators{
      MakeEigMat2d(footprints.obs_footprint_.length_ / 2.0, 0.0, 0.0,
                   footprints.obs_footprint_.width_ / 2.0)};
  const Eigen::Matrix<double, 2, 2> ego_generators_base{
      MakeEigMat2d(footprints.ego_footprint_.length_ / 2.0, 0.0, 0.0,
                   footprints.ego_footprint_.width_ / 2.0)};

  for (std::int64_t i = 0; i < ego_trajectory_global.cols(); ++i) {
    const Eigen::Vector3d& ego_traj_i{ego_trajectory_global.col(i)};
    const Eigen::Vector2d& obs_center_xy{obs_trajectory_global.col(i)};
    const Eigen::Vector2d ego_center_xy{ego_traj_i.x(), ego_traj_i.y()};
    const Eigen::Matrix2d ego_generators{
        Eigen::Rotation2Dd{ego_traj_i.z()}.matrix() * ego_generators_base};
    if (CheckIntersection(ego_center_xy, obs_center_xy, ego_generators,
                          obs_generators) ==
        CollisionStatus::kCollisionOccurred) {
      // fmt::print("[CO] Collision at ego ({}, {}) obs ({}, {})\n",
      // ego_center_xy.x(), ego_center_xy.y(), obs_center_xy.x(),
      // obs_center_xy.y());
      return CollisionStatus::kCollisionOccurred;
    }
    // fmt::print("[NO] Collision at ego ({}, {}) obs ({}, {})\n",
    // ego_center_xy.x(), ego_center_xy.y(), obs_center_xy.x(),
    // obs_center_xy.y());
  }
  return CollisionStatus::kNoCollisionOccurred;
}

struct MonteCarloConstraintOutput {
  bool risk_constraints_satisfied_;
  std::int64_t num_traj_samples_;
  std::int64_t num_collisions_;
};

struct FrsSampledInfo {
  FrsSelectInfo sel_inf_;
  double sampled_param_;
  std::int64_t num_traj_samples_;
  std::int64_t num_traj_collisions_;
  bool static_constraints_satisfied_;
  bool risk_constraints_satisfied_;
};

struct OptEConOut {
  std::int64_t num_traj_samples_;
  std::int64_t num_traj_collisions_;
  bool risk_constraints_satisfied_;
  double risk_;
  double risk_threshold_;
};

inline static OptEConOut
OptEConstraint(const Vehrs& param_sliced_frs,
               const std::vector<std::vector<double>>& all_mu_sigmas_local,
               const EnvFootprints& footprints,
               const std::int64_t num_trajectory_samples,
               const double risk_threshold, std::mt19937& gen) {
  const Eigen::Matrix2d obs_footprint_gens{
      MakeEigMat2d(footprints.obs_footprint_.length_, 0.0, 0.0,
                   footprints.obs_footprint_.width_)};
  const std::optional<double> buff_x{0.1};
  const std::optional<double> buff_y{0.1};
  // const std::optional<double> buff_x{std::nullopt};
  // const std::optional<double> buff_y{std::nullopt};

  double risk{0.0};
  std::int64_t total_check_iterations{0};
  std::int64_t total_num_collisions{0};

  for (const auto& mu_sigma_set : all_mu_sigmas_local) {
    if (risk > risk_threshold) {
      break;
    }
    if (mu_sigma_set.size() % 5 != 0) {
      throw std::runtime_error("Mu sigma set not divisible by 5");
    }
    if (mu_sigma_set.size() / 5 !=
        param_sliced_frs.zono_time_intervals_.size()) {
      throw std::runtime_error("mu sigma set has different num elts than FRS");
    }

    const std::int64_t num_mu_sigma_elts{
        static_cast<std::int64_t>(mu_sigma_set.size())};
    const std::int64_t num_zonos{static_cast<std::int64_t>(
        param_sliced_frs.zono_time_intervals_.size())};
    // Iterate over all mu sigmas
    for (int i = 0; i < num_mu_sigma_elts; i += 5) {
      if (risk > risk_threshold) {
        break;
      }
      // fmt::print("Obs Idx: {} of {}\n", i, num_mu_sigma_elts);
      //  Load mu sigma
      const double mu_x_local{mu_sigma_set.at(i + 0)};
      const double mu_y_local{mu_sigma_set.at(i + 1)};
      const double sigma_xx_local{mu_sigma_set.at(i + 2)};
      const double sigma_xy_local{mu_sigma_set.at(i + 3)};
      const double sigma_yy_local{mu_sigma_set.at(i + 4)};
      const IndividualMuSigma mu_sig_local{.mu_x_ = mu_x_local,
                                           .mu_y_ = mu_y_local,
                                           .sigma_1_ = sigma_xx_local,
                                           .sigma_2_ = sigma_xy_local,
                                           .sigma_4_ = sigma_yy_local};

      const int zono_idx = i / 5;
      // for (int zono_idx = 0; zono_idx < num_zonos; ++zono_idx) {
      //  fmt::print("Zono Idx: {} of {}\n", zono_idx, num_zonos);
      //  Compute ego center
      const double ego_c_x{param_sliced_frs.xy_centers_.at(zono_idx * 2 + 0)};
      const double ego_c_y{param_sliced_frs.xy_centers_.at(zono_idx * 2 + 1)};
      const Eigen::Vector2d ego_center{ego_c_x, ego_c_y};

      // Store ego generators
      Eigen::Matrix2Xd ego_gens;
      const auto num_ego_gens{param_sliced_frs.zono_sizes_.at(zono_idx)};
      ego_gens.resize(2, num_ego_gens);
      std::int64_t num_gen_offset{0};
      for (int prev_zono_idx = 0; prev_zono_idx + 1 < zono_idx;
           ++prev_zono_idx) {
        num_gen_offset += param_sliced_frs.zono_sizes_.at(prev_zono_idx);
      }

      for (int gen_idx = 0; gen_idx < num_ego_gens; ++gen_idx) {
        const std::int64_t gen_start{2 * (num_gen_offset + gen_idx)};
        const double gen_x{param_sliced_frs.zono_xy_.at(gen_start + 0)};
        const double gen_y{param_sliced_frs.zono_xy_.at(gen_start + 1)};
        const Eigen::Vector2d gen_xy{gen_x, gen_y};
        ego_gens.col(gen_idx) = gen_xy;
      }

      const auto t3{Tick()};
      // Sample and collision check
      std::int64_t num_collisions{0};
      for (int sample_no = 0; sample_no < num_trajectory_samples; ++sample_no) {
        const auto sampled_center{SampleMuSigma(mu_sig_local, gen)};
        ++total_check_iterations;
        if (CheckIntersectionArbitrary(ego_center, sampled_center, ego_gens,
                                       obs_footprint_gens, buff_x, buff_y) ==
            CollisionStatus::kCollisionOccurred) {
          ++num_collisions;
        }
      }
      total_num_collisions += num_collisions;
      const double risk_curr_pair{static_cast<double>(num_collisions) /
                                  static_cast<double>(num_trajectory_samples)};
      risk += risk_curr_pair;
      //}
    }
  }
  return {.num_traj_samples_ = total_check_iterations,
          .num_traj_collisions_ = total_num_collisions,
          .risk_constraints_satisfied_ = (risk <= risk_threshold),
          .risk_ = risk,
          .risk_threshold_ = risk_threshold};
}

inline static double SampleParamVal(const Vehrs& frs_individ,
                                    std::mt19937& gen) {
  const auto min_traj_param{frs_individ.GetTrajParamMin()};
  const auto max_traj_param{frs_individ.GetTrajParamMax()};
  return SampleUniform(min_traj_param, max_traj_param, gen);
}

inline static std::vector<RiskRtdPlanningOutputs>
MonteCarloDiscreteSamplingOptE(const ::roahm::FrsTotal& frs,
                               const RiskProblemDescription& risk_inputs,
                               const std::uint32_t randomization_seed,
                               const std::int64_t num_trajectory_samples,
                               const ::roahm::TimeMicroseconds dt,
                               const EnvFootprints footprints,
                               const double risk_threshold) {
  const int param_samples_per_traj{1};
  fmt::print("[BEGIN] MonteCarloDiscreteSampling Ego: {} x {} (Lat x Long)\n",
             footprints.obs_footprint_.length_,
             footprints.ego_footprint_.length_);

  const auto& rover_state{risk_inputs.rover_state_};

  std::mt19937 gen{randomization_seed};

  const auto all_valid_frs_infos{FindAllValidFrses(frs, rover_state)};
  fmt::print("Number of valid FRSES: {} from FRS with {} megas\n",
             all_valid_frs_infos.size(), frs.megas_.size());

  std::vector<RiskRtdPlanningOutputs> out;
  std::vector<FrsSampledInfo> frs_sampled_info{};
  for (const auto& valid_frs_info : all_valid_frs_infos) {
    // Load the FRS
    const auto& frs_individ{frs.GetVehrs(valid_frs_info)};

    for (int param_sample_idx = 0; param_sample_idx < param_samples_per_traj;
         ++param_sample_idx) {
      // Sample the trajectory parameter from the uniform distribution [min,
      // max]
      const double sampled_param_val{SampleParamVal(frs_individ, gen)};

      const auto mirror{valid_frs_info.mirror_};

      // Load initial state
      const double u0{rover_state.u_};
      const double v0{rover_state.v_};
      const double r0{rover_state.r_};

      // Slice the FRS
      const auto state_sliced_frs{frs_individ.SliceAt(u0, v0, r0)};
      const auto param_sliced_frs{
          frs_individ.SliceAtParam(u0, v0, r0, sampled_param_val)};

      auto always_non_risky{risk_inputs.always_non_risky_obs_.RelativeTo(
          rover_state.GetXYH(), mirror)};
      auto always_risky{risk_inputs.always_risky_obs_};
      for (const auto& maybe_risky_obs : risk_inputs.maybe_risky_obs_) {
        if (maybe_risky_obs.dyn_obs_.IsStaticBoundary()) {
          always_non_risky.PushObs(maybe_risky_obs.dyn_obs_.RelativeTo(
              rover_state.GetXYH(), mirror));
        } else {
          always_risky.push_back(maybe_risky_obs.mu_sigmas_);
        }
      }

      for (int i = 0; i < always_non_risky.GetNumObs(); ++i) {
        auto curr_obs{always_non_risky.GetSingleDynObs(i)};
        fmt::print("OBS {} ({}, {}), LW ({}, {}), V({})\n", i,
                   curr_obs.CenterX0(), curr_obs.CenterY0(), curr_obs.Length(),
                   curr_obs.Width(), curr_obs.Velocity());
      }

      // Generate local frame waypoints
      const auto waypoint_local_mirror_accounted{
          risk_inputs.waypoint_global_no_mirror_
              .RelativeToEgoFrameNoMirror(risk_inputs.rover_state_.GetXYH())
              .TakeMirrorIntoAccount(mirror)};
      const auto waypoint_heuristic_adjusted{
          waypoint_local_mirror_accounted.HeuristicAdjust(
              frs_individ.GetManuType(), rover_state.u_)};
      const auto waypoint_cost_fcn = ::roahm::WaypointCostFromVehrsAndWp(
          state_sliced_frs, waypoint_heuristic_adjusted);
      const auto final_location =
          waypoint_cost_fcn.GetPointAt(sampled_param_val);
      const double cost{waypoint_cost_fcn.EvaluateAt(sampled_param_val).cost_};

      // Check static constraints
      auto static_constraints{::roahm::fl_zono_constraint::GetFlZonoConstraint(
          state_sliced_frs, always_non_risky)};
      const bool static_satisfied{
          static_constraints.ConstraintsSatisfied(sampled_param_val)};
      if (not static_satisfied) {
        fmt::print("STATIC NOT SATISFIED\n");
        for (int i = 0; i < always_non_risky.GetNumObs(); ++i) {
          const auto curr_obs{always_non_risky.GetSingleDynObs(i)};
          ::roahm::FlZonoObsSet obs_set_curr;
          obs_set_curr.PushObs(curr_obs);
          auto curr_cons{fl_zono_constraint::GetFlZonoConstraint(
              state_sliced_frs, obs_set_curr)};
          if (curr_cons.ConstraintsSatisfied(sampled_param_val)) {
            fmt::print("CONSTRAINT {}: SATISFIED\n", i);
          } else {
            fmt::print("CONSTRAINT {}: NOT SATISFIED\n", i);
          }
        }
        frs_sampled_info.push_back({.sel_inf_ = valid_frs_info,
                                    .sampled_param_ = sampled_param_val,
                                    .num_traj_samples_ = 0,
                                    .num_traj_collisions_ = 0,
                                    .static_constraints_satisfied_ = false,
                                    .risk_constraints_satisfied_ = false});
      } else {
        fmt::print("STATIC SATISFIED\n");
        std::vector<std::vector<double>> all_mu_sigmas;
        for (const auto& always_risky_mu_sigma_multi : always_risky) {
          all_mu_sigmas.push_back(frs_individ.ExtractMuSigma(
              always_risky_mu_sigma_multi, rover_state.GetXYH(), mirror));
        }
        fmt::print("NUM MU SIGMAS: {}\n", all_mu_sigmas.size());
        const auto opt_e_ret{OptEConstraint(param_sliced_frs, all_mu_sigmas,
                                            footprints, num_trajectory_samples,
                                            risk_threshold, gen)};
        fmt::print("RISK: {}\n", opt_e_ret.risk_);
        frs_sampled_info.push_back(
            {.sel_inf_ = valid_frs_info,
             .sampled_param_ = sampled_param_val,
             .num_traj_samples_ = opt_e_ret.num_traj_samples_,
             .num_traj_collisions_ = opt_e_ret.num_traj_collisions_,
             .static_constraints_satisfied_ = true,
             .risk_constraints_satisfied_ =
                 opt_e_ret.risk_constraints_satisfied_});
      }

      const auto last_sampled{frs_sampled_info.back()};
      const bool last_successful{last_sampled.static_constraints_satisfied_ and
                                 last_sampled.risk_constraints_satisfied_};
      fmt::print("OPT E: {}/{} col\n", last_sampled.num_traj_collisions_,
                 last_sampled.num_traj_samples_);
      const RiskRtdPlanningOutputs out_curr{
          .ran_ = true,
          .solve_succeeded_ = last_successful,
          .found_feasible_ = last_successful,
          .sel_inf_ = valid_frs_info,
          .final_param_ = sampled_param_val,
          .final_cost_ = cost,
          .u0_ = rover_state.u_,
          .final_location_ = final_location,
          .waypoint_with_heuristic_ = waypoint_heuristic_adjusted,
          .waypoint_local_frame_ = waypoint_local_mirror_accounted,
      };
      out.push_back(out_curr);
    }
  }

  return out;
}

inline static
OptEConOut MonteCarloConstraint27(
               const Vehrs& param_sliced_frs,
               const bool mirror,
               const std::vector<MuSigmaMulti>& always_risky,
               const EnvFootprints& footprints,
               const std::int64_t num_trajectory_samples,
               const TimeMicroseconds& dt,
               const RoverState& z0,
               const double risk_threshold, 
               std::mt19937& gen) {
  // Slice the FRS at the chosen parameter
  const ::roahm::TimeMicroseconds trajectory_duration{
      param_sliced_frs.TrajectoryDuration()};

  const auto ego_traj_global{
      GetEgoTrajectoryGlobalFrame(param_sliced_frs, dt, z0.GetXYH())};
  // Sample each obstacle's trajectory and collision check
  std::int64_t num_collisions{0};
  std::vector<std::vector<IndividualMuSigma>> all_individ_mu_sigmas;
  for (std::int64_t i = 0; i < num_trajectory_samples; ++i) {
    const double curr_risk_lb{static_cast<double>(num_collisions) /
                              static_cast<double>(num_trajectory_samples)};
    if (curr_risk_lb > risk_threshold) {
      break;
    }
    for (const auto& mu_sigmas : always_risky) {
      const auto sampled_trajectory{SampleSingleTrajectory(
          mu_sigmas, trajectory_duration, dt, gen, z0.GetXYH(), mirror)};
      const auto cc_ret{CollisionCheckTrajectory(
          ego_traj_global, sampled_trajectory, footprints)};
      if (cc_ret == CollisionStatus::kCollisionOccurred) {
        ++num_collisions;
        break;
      }
    }
  }
  const double risk{static_cast<double>(num_collisions) /
                    static_cast<double>(num_trajectory_samples)};
  const bool risk_constraints_satisfied{risk < risk_threshold};
  return {.num_traj_samples_ = num_trajectory_samples,
          .num_traj_collisions_ = num_collisions,
          .risk_constraints_satisfied_ = risk_constraints_satisfied,
          .risk_ = risk,
          .risk_threshold_ = risk_threshold};
}

inline static std::vector<RiskRtdPlanningOutputs> MonteCarloDiscreteSampling27(
    const ::roahm::FrsTotal& frs, 
    const RiskProblemDescription& risk_inputs,
    const std::uint32_t randomization_seed,
    const std::int64_t num_trajectory_samples,
    const ::roahm::TimeMicroseconds dt, 
    const EnvFootprints footprints,
    const double risk_threshold,
    const int num_param_samples) {
  fmt::print("[BEGIN] MonteCarloDiscreteSampling Ego: {} x {} (Lat x Long)\n",
             footprints.obs_footprint_.length_,
             footprints.ego_footprint_.length_);
  const auto waypoint_local_no_mirror{
      risk_inputs.waypoint_global_no_mirror_.RelativeToEgoFrameNoMirror(
          risk_inputs.rover_state_.GetXYH())};

  const auto& rover_state{risk_inputs.rover_state_};

  std::mt19937 gen{randomization_seed};

  const auto all_valid_frs_infos{FindAllValidFrses(frs, rover_state)};
  fmt::print("Number of valid FRSES: {} from FRS with {} megas\n",
             all_valid_frs_infos.size(), frs.megas_.size());

  std::vector<RiskRtdPlanningOutputs> out;
  std::vector<FrsSampledInfo> frs_sampled_info{};
  for (const auto& valid_frs_info : all_valid_frs_infos) {
    // Load the FRS
    const auto& frs_individ{frs.GetVehrs(valid_frs_info)};
    // Sample the trajectory parameter from the uniform distribution [min, max]
    const auto min_traj_param{frs_individ.GetTrajParamMin()};
    const auto max_traj_param{frs_individ.GetTrajParamMax()};
    const auto mirror{valid_frs_info.mirror_};
    for (int param_sample_no = 0; param_sample_no < num_param_samples; ++param_sample_no) {
      const double sampled_param_val{
          SampleUniform(min_traj_param, max_traj_param, gen)};

      // Construct Obs
      ::roahm::FlZonoObsSet always_non_risky{risk_inputs.always_non_risky_obs_.RelativeTo(rover_state.GetXYH(), mirror)};
      auto always_risky{risk_inputs.always_risky_obs_};
      for (const auto& maybe_risky_obs : risk_inputs.maybe_risky_obs_) {
        const auto curr_dyn_obs{maybe_risky_obs.dyn_obs_};
        if (curr_dyn_obs.IsStaticBoundary()) {
          always_non_risky.PushObs(curr_dyn_obs.RelativeTo(rover_state.GetXYH(), mirror));
        } else {
          always_risky.push_back(maybe_risky_obs.mu_sigmas_);
        }
      }


      // Load initial state
      const double u0{rover_state.u_};
      const double v0{rover_state.v_};
      const double r0{rover_state.r_};

      // Slice the FRS
      const auto state_sliced_frs{frs_individ.SliceAt(u0, v0, r0)};
      const auto param_sliced_frs{
          frs_individ.SliceAtParam(u0, v0, r0, sampled_param_val)};

      // Generate local frame waypoints
      const auto waypoint_local_mirror_accounted{
          waypoint_local_no_mirror.TakeMirrorIntoAccount(mirror)};
      const auto waypoint_heuristic_adjusted{
          waypoint_local_mirror_accounted.HeuristicAdjust(
              frs_individ.GetManuType(), rover_state.u_)};
      const auto waypoint_cost_fcn = ::roahm::WaypointCostFromVehrsAndWp(
          state_sliced_frs, waypoint_heuristic_adjusted);
      const auto final_location = waypoint_cost_fcn.GetPointAt(sampled_param_val);
      const double cost{waypoint_cost_fcn.EvaluateAt(sampled_param_val).cost_};

      // Check static constraints
      auto static_constraints{::roahm::fl_zono_constraint::GetFlZonoConstraint(
          state_sliced_frs, always_non_risky)};
      const bool static_satisfied{
          static_constraints.ConstraintsSatisfied(sampled_param_val)};
      if (not static_satisfied) {
        frs_sampled_info.push_back({.sel_inf_ = valid_frs_info,
                                    .sampled_param_ = sampled_param_val,
                                    .num_traj_samples_ = 0,
                                    .num_traj_collisions_ = 0,
                                    .static_constraints_satisfied_ = false,
                                    .risk_constraints_satisfied_ = false});
      } else {
        const auto constraint_eval{MonteCarloConstraint27(
            param_sliced_frs, mirror, always_risky, footprints,
            num_trajectory_samples, dt, rover_state, risk_threshold, gen)};
        frs_sampled_info.push_back(
            {.sel_inf_ = valid_frs_info,
             .sampled_param_ = sampled_param_val,
             .num_traj_samples_ = constraint_eval.num_traj_samples_,
             .num_traj_collisions_ = constraint_eval.num_traj_samples_,
             .static_constraints_satisfied_ = true,
             .risk_constraints_satisfied_ =
                 constraint_eval.risk_constraints_satisfied_});
      }

      const auto last_sampled{frs_sampled_info.back()};
      const bool last_successful{last_sampled.static_constraints_satisfied_ and
                                 last_sampled.risk_constraints_satisfied_};
      fmt::print("{}/{} col\n", last_sampled.num_traj_collisions_,
                 last_sampled.num_traj_samples_);
      const RiskRtdPlanningOutputs out_curr{
          .ran_ = true,
          .solve_succeeded_ = last_successful,
          .found_feasible_ = last_successful,
          .sel_inf_ = valid_frs_info,
          .final_param_ = sampled_param_val,
          .final_cost_ = cost,
          .u0_ = rover_state.u_,
          .final_location_ = final_location,
          .waypoint_with_heuristic_ = waypoint_heuristic_adjusted,
          .waypoint_local_frame_ = waypoint_local_mirror_accounted,
      };
      out.push_back(out_curr);
    }
  }

  return out;
}

} // namespace roahm::risk_comparisons
#endif // ROAHM_MONTE_CARLO_OPT_E_HPP_
