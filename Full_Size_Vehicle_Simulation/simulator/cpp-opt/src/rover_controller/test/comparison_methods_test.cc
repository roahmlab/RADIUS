#include "comparison_methods.hpp"
#include "dyn_obs.hpp"
#include "fl_zono_obs_set.hpp"
#include "mu_sigma_multi.hpp"
#include "risk_rtd.hpp"
#include "risk_rtd_ipopt_problem.hpp"
#include "ros/time.h"
#include "waypoint.hpp"
#include <gtest/gtest.h>
#include <random>

TEST(AblationCompile, Ablate) {
  const std::string frs_fname{
      "/data/cpp_processed_CUDA_FRS_30-Aug-2022_lessFRS_3rnd_grid24_t_fix.DBG"};
  const auto frs_full{roahm::LoadFrsBinary(frs_fname)};
  const ::roahm::RoverState init_state{0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 20.0};
  std::array<std::vector<double>, 16> mu_sigmas_dat;
  for (int i = 0; i < mu_sigmas_dat.size(); ++i) {
    for (int j = 0; j < static_cast<std::int64_t>(13.0 / 0.01); ++j) {
      const double mu_x{0.0 + 0.01 * j * 20.0};
      const double mu_y{2.0};
      const double sigma_xx{1.0};
      const double sigma_xy{0.5};
      const double sigma_yy{0.3};
      mu_sigmas_dat.at(i).push_back(mu_x);
      mu_sigmas_dat.at(i).push_back(mu_y);
      mu_sigmas_dat.at(i).push_back(sigma_xx);
      mu_sigmas_dat.at(i).push_back(sigma_xy);
      mu_sigmas_dat.at(i).push_back(sigma_yy);
    }
  }

  const roahm::TimeMicroseconds dt{roahm::TimeMicroseconds::Get10ms()};
  roahm::MuSigmaMulti mu_sigma_mult{mu_sigmas_dat, dt};

  ::roahm::FlZonoObsSet always_non_risky;
  always_non_risky.PushObs(
      ::roahm::DynObs(0.0, 3.5, 0.0, 20.0, 1.0, 1.0,
                      ::roahm::DynObs::ObstacleType::kStaticBoundary));

  const ::roahm::RiskProblemDescription problem_description{
      .rover_state_{init_state},
      .maybe_risky_obs_{},
      .always_risky_obs_{mu_sigma_mult},
      .always_non_risky_obs_{always_non_risky},
      .waypoint_global_no_mirror_{}};
  //::roahm::MonteCarloDiscreteSampling(frs_full, problem_description, 0);

  const double mu_x{-18};
  const double mu_y{-5.0};
  const double sigma_xx{0.001};
  const double sigma_xy{-0.02};
  const double sigma_yy{30.0};
  std::random_device rd;
  std::mt19937 gen(rd());
  for (int i = 0; i < 1000; ++i) {
    const auto ret{::roahm::risk_comparisons::GetSample(
        mu_x, mu_y, sigma_xx, sigma_xy, sigma_yy, gen)};
    fmt::print("{},{}\n", ret.x(), ret.y());
  }
  const std::uint32_t rand_seed{0};
  const std::int64_t num_trajectory_samples{1000};

  const roahm::WaypointGlobalNoMirror wp_global{{100.0, 0.0, 0.0}};
  const roahm::risk_comparisons::EnvFootprints footprints{
      ::roahm::risk_comparisons::EnvFootprints::Default()};
  /*
  roahm::risk_comparisons::MonteCarloDiscreteSampling(
      frs_full, problem_description, rand_seed, num_trajectory_samples,
      ::roahm::TimeMicroseconds::Get10ms(), footprints, wp_global, 0.05);

  auto tmp_apps{roahm::risk_rtd_ipopt_problem::GetRiskIpoptApps(1)};
  */

  roahm::risk_comparisons::MonteCarloDiscreteSamplingOptE(
      frs_full, problem_description, rand_seed, num_trajectory_samples, dt,
      footprints, 0.05);

  ::roahm::RiskRtd risk_rtd{frs_fname};
  risk_rtd.RunPlanningIteration<
      ::roahm::risk_rtd_ipopt_problem::RiskRtdIpoptProblem>(
      problem_description, true, 0.05);
}
