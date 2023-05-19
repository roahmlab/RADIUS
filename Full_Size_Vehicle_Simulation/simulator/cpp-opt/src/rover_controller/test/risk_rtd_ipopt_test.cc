#include <gtest/gtest.h>

#include <fstream>
#include <string>

#include "frs_io.hpp"
#include "frs_loader.hpp"
#include "ipopt_string_utils.hpp"
#include "risk_rtd_ipopt_problem.hpp"
#include "timing_util.hpp"

std::vector<double> BuildMuSigma(const int num_zonos) {
  std::vector<double> mu_sigma;
  const double mu_1 = 0.84;
  const double mu_2 = 0.91;
  const double sigma_1 = 1.42;
  const double sigma_2 = 0.21;
  const double sigma_4 = 0.31;
  for (int i = 0; i < num_zonos; ++i) {
    mu_sigma.push_back(mu_1);
    mu_sigma.push_back(mu_2);
    mu_sigma.push_back(sigma_1);
    mu_sigma.push_back(sigma_2);
    mu_sigma.push_back(sigma_4);
  }
  return mu_sigma;
}

TEST(RiskRtdIpoptTest, Test) {
  using ::roahm::GetDeltaS;
  using ::roahm::Tick;
  const std::string frs_in_fname{"/data/frs_bin_processed.txt"};
  std::ifstream file_in(frs_in_fname, std::ios::binary);
  const ::roahm::FrsTotal frs_in{
      ::roahm::ReadFromBinFile<::roahm::FrsTotal>(file_in)};
  file_in.close();

  const auto ipopt_apps = ::roahm::risk_rtd_ipopt_problem::GetRiskIpoptApps(1);

  auto& app = ipopt_apps.at(0);
  if (app->Initialize() != Ipopt::Solve_Succeeded) {
    std::cout << "\n\n*** Error during IPOPT initialization!" << std::endl;
  }
  const auto vehrs = frs_in.megas_.at(0).au_.at(0);
  const std::array<double, 3> u0v0r0_slice_beta{0.0, 0.0, 0.0};
  const auto mu_sigma = BuildMuSigma(vehrs.cuda_info_.num_zono_);
  const std::vector<std::vector<double>> all_mu_sigmas{mu_sigma};
  ::roahm::OptimizationInputs opt_inputs;
  opt_inputs.u0_meters_per_second_ = 0.0;
  opt_inputs.mu_sigmas_ = all_mu_sigmas;
  opt_inputs.use_fixed_risk_threshold_ = true;
  opt_inputs.risk_threshold_ = 0.05;
  const auto t0 = Tick();
  Ipopt::SmartPtr<Ipopt::TNLP> mynlp =
      new ::roahm::risk_rtd_ipopt_problem::RiskRtdIpoptProblem(
          vehrs, u0v0r0_slice_beta, opt_inputs);
  const auto t1 = Tick();
  const Ipopt::ApplicationReturnStatus status = app->OptimizeTNLP(mynlp);
  const auto t2 = Tick();
  std::cout << "Setup Took: " << GetDeltaS(t1, t0) << std::endl;
  std::cout << "Optimizing Took: " << GetDeltaS(t2, t1) << std::endl;
  std::cout << "Result: " << ::roahm::ToString(status) << std::endl;
  const bool ipopt_success = status == Ipopt::Solve_Succeeded or
                             status == Ipopt::Solved_To_Acceptable_Level;

  EXPECT_EQ(1, 1);
}
