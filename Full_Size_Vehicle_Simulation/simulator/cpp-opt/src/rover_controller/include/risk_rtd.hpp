#ifndef ROAHM_RISK_RTD_HPP_
#define ROAHM_RISK_RTD_HPP_

#include <chrono>
#include <limits>
#include <map>
#include <optional>
#include <thread>
#include <type_traits>
#include <unistd.h>
#include <vector>

#include "dyn_obs.hpp"
#include "fl_zono_obs_set.hpp"
#include "frs_loader.hpp"
#include "frs_select_info.hpp"
#include "frs_total.hpp"
#include "ipopt_string_utils.hpp"
#include "manu_type.hpp"
#include "monte_carlo_ipopt.hpp"
#include "monte_carlo_27_ipopt.hpp"
#include "mu_sigma_multi.hpp"
#include "risk_problem_description.hpp"
#include "risk_rtd_ipopt_problem.hpp"
#include "rover_state.hpp"
#include "waypoint.hpp"

namespace roahm {

/// Contains both the mu sigma and deterministic representation of an obstacle,
/// whether to treat it as a risk constraint or a non-risk constraint is
/// determined at runtime.
struct RiskTrajectoryInfo {
  FrsSelectInfo frs_info_;
  double param_val_;
};

struct RiskRtdPlanningOutputsIpopt : public RiskRtdPlanningOutputs {
  Ipopt::ApplicationReturnStatus ipopt_status_;
  std::optional<Ipopt::SolverReturn> final_solver_status_;
};

class RiskRtd {
private:
  ::roahm::FrsTotal frs_;
  RiskRtd() = delete;

public:
  RiskRtd(const std::string& binary_file_name)
      : frs_{::roahm::LoadFrsBinary(binary_file_name)} {}

  template <typename OptProblemType>
  std::vector<RiskRtdPlanningOutputsIpopt>
  RunPlanningIteration(const RiskProblemDescription& risk_inputs, const bool use_fixed_risk_threshold, const double risk_threshold) const;

  FrsTotal GetFrsTotal() const;
};

template <typename... Args>
void PlanDbg(const fmt::format_string<Args...>&& format_string,
             Args&&... args) {
  fmt::print("[PLAN IT] {}",
             fmt::format(format_string, std::forward<Args>(args)...));
}

template <typename OptProblemType>
std::vector<RiskRtdPlanningOutputsIpopt>
RiskRtd::RunPlanningIteration(const RiskProblemDescription& risk_inputs, const bool use_fixed_risk_threshold, const double risk_threshold) const {
  const auto iter_t0 = Tick();
  PlanDbg("Running planning iteration begin\n");
  PlanDbg("Finding valid FRSes\n");
  const auto sel = FindAllValidFrses(frs_, risk_inputs.rover_state_);
  PlanDbg("Found {} valid FRSes out of {}\n", sel.size(), frs_.GetNumFrses());

  const int num_apps = 6;
  const int num_problems = sel.size();
  struct SingleOptRet {
    bool feasible_;
    ManuType manu_type_;
    double param_val_;
  };
  PlanDbg("Setting up {} risk apps\n", num_apps);
  auto apps = ::roahm::risk_rtd_ipopt_problem::GetRiskIpoptApps(num_apps);
  PlanDbg("Setting opt_rets\n");
  std::vector<SingleOptRet> opt_rets(num_problems);
  PlanDbg("Setting atomics\n");
  std::atomic<int> iter_ct{0};
  std::atomic<int> success_ct{0};
  PlanDbg("Setting nlps\n");
  std::vector<Ipopt::SmartPtr<Ipopt::TNLP>> nlps;
  PlanDbg("Reserving nlps\n");
  nlps.reserve(num_problems);
  for (int problem_no = 0; problem_no < num_problems; ++problem_no) {
    const auto& info = sel.at(problem_no);
    const auto vehrs = frs_.GetVehrs(info).SliceAt(risk_inputs.rover_state_.u_,
                                                   risk_inputs.rover_state_.v_,
                                                   risk_inputs.rover_state_.r_);
    const auto opt_inputs = GetOptInputs(frs_, info, risk_inputs, use_fixed_risk_threshold, risk_threshold);
    const std::array<double, 3> u0v0r0_slice_beta =
        vehrs.GetU0V0R0SliceBeta(risk_inputs.rover_state_);
    if constexpr (std::is_same_v<
                      OptProblemType,
                      ::roahm::risk_rtd_ipopt_problem::RiskRtdIpoptProblem>) {
      nlps.emplace_back(
          new ::roahm::risk_rtd_ipopt_problem::RiskRtdIpoptProblem{
              vehrs, u0v0r0_slice_beta, opt_inputs});
    } else if (std::is_same_v<
                   OptProblemType,
                   ::roahm::monte_carlo_ipopt::MonteCarloIpoptProblem>) {
      const int num_traj_samples{100};
      const std::uint32_t rand_seed{0};
      nlps.emplace_back(new ::roahm::monte_carlo_ipopt::MonteCarloIpoptProblem{
          vehrs, u0v0r0_slice_beta, opt_inputs, rand_seed, frs_.GetVehrs(info),
          risk_inputs.rover_state_, risk_comparisons::EnvFootprints::Default(),
          num_traj_samples});
    } else if (std::is_same_v<OptProblemType,
    ::roahm::monte_carlo_27_ipopt::MonteCarlo27IpoptProblem>) {
      const int num_traj_samples{100};
      const std::uint32_t rand_seed{0};
      nlps.emplace_back(new ::roahm::monte_carlo_27_ipopt::MonteCarlo27IpoptProblem{
        vehrs, u0v0r0_slice_beta, opt_inputs, rand_seed, frs_.GetVehrs(info),
        risk_inputs.rover_state_, risk_comparisons::EnvFootprints::Default(),
        num_traj_samples,
        info.mirror_,
        risk_inputs
      });
    }
    // static_assert(kIsRiskRtd);
  }

  PlanDbg("Setting up problem function\n");
  PlanDbg("Num Apps (Threads): {}\n", num_apps);
  PlanDbg("Num Problems: {}\n", num_problems);
  std::vector<RiskRtdPlanningOutputsIpopt> full_opt_info_all_problems_;
  full_opt_info_all_problems_.resize(num_problems);
  const double u0 = risk_inputs.rover_state_.u_;
  auto run_problems = [&apps, &iter_ct, &nlps, &success_ct, &sel,
                       &full_opt_info_all_problems_, &u0,
                       num_problems](const int app_idx) -> void {
    auto& ipopt_app = apps.at(app_idx);
    while (true) {
      const auto t0 = Tick();
      const auto problem_no = iter_ct++;
      if (problem_no >= num_problems) {
        break;
      }
      PlanDbg("[START] Problem {} on thread {}\n", problem_no, app_idx);
      if (ipopt_app->Initialize() != Ipopt::Solve_Succeeded) {
        fmt::print("\n\n***** Error during IPOPT initialization!\n");
      }
      auto curr_nlp_base = nlps.at(problem_no);
      const auto t1 = Tick();
      const Ipopt::ApplicationReturnStatus status =
          ipopt_app->OptimizeTNLP(curr_nlp_base);
      const auto t2 = Tick();
      const auto* curr_nlp =
          dynamic_cast<OptProblemType*>(GetRawPtr(curr_nlp_base));
      std::optional<double> final_cost{std::nullopt};
      std::optional<double> final_param{std::nullopt};
      bool found_feasible{false};
      bool solve_succeeded{false};
      // TODO NOTE: currently == solve succeeded
      const auto final_solver_status = curr_nlp->GetSolverReturnStatus();
      found_feasible = curr_nlp->FoundFeasible();
      final_param = curr_nlp->GetFeasibleParam();
      final_cost = curr_nlp->GetFeasibleCost();
      const std::optional<PointXYH> final_location{
          curr_nlp->GetFeasibleLocation()};
      const Waypoint<true, true> waypoint_with_heuristic{
          curr_nlp->GetWaypointHeuristicAdjusted()};
      const Waypoint<true, true> waypoint_local_mirror_accounted{
          curr_nlp->GetWaypointLocalMirrorAccounted()};
      // TODO or: status == Ipopt::ApplicationReturnStatus::Solve_Succeeded or
      // ??
      if (found_feasible) {
        // if (status !=
        // Ipopt::ApplicationReturnStatus::Infeasible_Problem_Detected) {
        ++success_ct;
      }
      RiskRtdPlanningOutputsIpopt full_opt_inf{
          RiskRtdPlanningOutputs{
              .ran_ = true,
              .solve_succeeded_ = solve_succeeded,
              .found_feasible_ = found_feasible,
              .sel_inf_ = sel.at(problem_no),
              .final_param_ = final_param,
              .final_cost_ = final_cost,
              .u0_ = u0,
              .final_location_ = final_location,
              .waypoint_with_heuristic_ = waypoint_with_heuristic,
              .waypoint_local_frame_ = waypoint_local_mirror_accounted,
          },
          status,
          final_solver_status,
      };
      full_opt_info_all_problems_.at(problem_no) = full_opt_inf;
      const std::string solver_ret_status =
          final_solver_status.has_value()
              ? ToString(final_solver_status.value())
              : "nullopt";

      PlanDbg("[ END ] Problem {} on thread {} with status {} and feasibility "
              "{} [t1t0: {}] "
              "[t2t1: {}]\n",
              problem_no, app_idx, ToString(status), found_feasible,
              GetDeltaS(t1, t0), GetDeltaS(t2, t1));
    }
  };

  const auto iter_t1 = Tick();
  if (num_apps > 1) {
    PlanDbg("Creating threads...\n");
    std::vector<std::thread> threads;
    for (int i = 0; i < num_apps; ++i) {
      PlanDbg("Emplacing {}\n", i);
      threads.emplace_back(std::thread(run_problems, i));
    }

    PlanDbg("Joining all...\n");
    for (auto& problem_thread : threads) {
      PlanDbg("Joining...\n");
      problem_thread.join();
    }
  } else {
    run_problems(0);
  }
  const auto iter_t2 = Tick();
  PlanDbg("Joined all threads\n");
  // run_problems(0);

  // for (int i = 0; i < num_problems; ++i) {
  //   const auto& info = sel.at(iter_ct++);
  //   const auto vehrs =
  //   frs_.GetVehrs(info).SliceAt(risk_inputs.rover_state_.u_,
  //                                                  risk_inputs.rover_state_.v_,
  //                                                  risk_inputs.rover_state_.r_);
  //   // TODO should this be clamped?
  //   const auto opt_inputs = GetOptInputs(frs_, info, risk_inputs);
  //   const std::array<double, 3> u0v0r0_slice_beta =
  // vehrs.GetU0V0R0SliceBeta(risk_inputs.rover_state_); auto& app =
  // apps.front();
  //   if (app->Initialize() != Ipopt::Solve_Succeeded) {
  //     std::cout << "\n\n*** Error during IPOPT initialization!" << std::endl;
  //   }
  //   Ipopt::SmartPtr<Ipopt::TNLP> mynlp =
  //       new ::roahm::risk_rtd_ipopt_problem::RiskRtdIpoptProblem{
  //           vehrs, u0v0r0_slice_beta, opt_inputs};
  //   const Ipopt::ApplicationReturnStatus status = app->OptimizeTNLP(mynlp);
  // }

  // auto t3 = Tick();

  // std::cout << "Finding took:    " << GetDeltaS(t2, t1) << "[s]" <<
  // std::endl; std::cout << "Waypoint [Local, No Mirror]: "
  //           << risk_inputs.waypoint_local_no_mirror_ << std::endl;
  // std::cout << "Mu Sigma took:   " << GetDeltaS(t3, t2) << "[s]" <<
  // std::endl; std::cout << "Rover state: " << risk_inputs.rover_state_.x_ <<
  // std::endl; std::cout << "U0 Idx: "
  //           << frs_.SelectU0IdxNoClamp(risk_inputs.rover_state_.u_)
  //           << std::endl;
  // std::cout << "Found " << sel.size() << " valid FRSes" << std::endl;

  // for (int i = 0; i < 10; ++i) {
  //   std::cout << "MuSig[" << i << "]: "
  //             << risk_inputs.maybe_risky_obs_.back()
  //                    .mu_sigmas_.mu_sigma_10ms_full_.at(i)
  //             << std::endl;
  //   std::cout << "MuSig[" << i << "]: "
  //             << risk_inputs.maybe_risky_obs_.back()
  //                    .mu_sigmas_.mu_sigma_20ms_full_.at(i)
  //             << std::endl;
  // }
  // RiskRtdPlanningOutputs outputs;
  // outputs.total_time_seconds_ = GetDeltaS(t3, t1);
  // return outputs;
  const auto iter_t3 = Tick();
  PlanDbg("Generating planning iteration output\n");
  PlanDbg("Planning iteration returning \n");
  PlanDbg("Total Plan Time: {}s\n", GetDeltaS(iter_t3, iter_t0));
  PlanDbg("T1-T0: {}s\n", GetDeltaS(iter_t1, iter_t0));
  PlanDbg("T2-T1: {}s\n", GetDeltaS(iter_t2, iter_t1));
  PlanDbg("T3-T2: {}s\n", GetDeltaS(iter_t3, iter_t2));
  return full_opt_info_all_problems_;
}

} // namespace roahm

#endif // ROAHM_RISK_RTD_HPP_
