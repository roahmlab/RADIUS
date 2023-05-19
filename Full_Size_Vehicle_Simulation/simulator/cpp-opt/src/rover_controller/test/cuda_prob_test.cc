#include <gtest/gtest.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include "frs_io.hpp"
#include "frs_loader.hpp"
#include "pre_slice.hpp"
#include "prob_integ.hpp"
#include "risk_constraint.hpp"
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

TEST(FakeTest, TestFake) {
  const std::string frs_in_fname{"/data/frs_bin_processed.txt"};
  std::ifstream file_in(frs_in_fname, std::ios::binary);
  const ::roahm::FrsTotal frs_in{
      ::roahm::ReadFromBinFile<::roahm::FrsTotal>(file_in)};
  file_in.close();

  const auto& cuda_info = frs_in.megas_.at(0).au_.at(0).cuda_info_;
  const std::array<double, 3> u0v0r0_slice_beta{0.13, 0.77, 0.04};
  const std::vector<double> mu_sigma = BuildMuSigma(cuda_info.num_zono_);
  const auto outputs = ::roahm::pre_slice::PreSlice(
      cuda_info, u0v0r0_slice_beta, mu_sigma.data());
  double p = 0.01;
  p = (p - cuda_info.cg_p_.at(0)) / cuda_info.cg_p_.at(1);
  std::cout << "Num Zonos: " << cuda_info.num_zono_ << std::endl;
  std::cout << "P: " << p << std::endl;
  const ::roahm::ProbIntegrationInputs prob_inputs{
      cuda_info.num_zono_,
      cuda_info.grid_x0_.data(),
      cuda_info.grid_y0_.data(),
      cuda_info.grid_dx_.data(),
      cuda_info.grid_dy_.data(),
      cuda_info.block_inzono_list_.data(),
      cuda_info.rot_angle_.data(),
      mu_sigma.data(),
      cuda_info.cg_p_.data(),
      cuda_info.g_p_x_.data(),
      cuda_info.grid_size_,
      p,
      outputs.H1_.get(),
      outputs.H2_.get(),
      outputs.H4_.get()};
  const auto prob_outputs = ::roahm::prob_integ::ProbIntegration(prob_inputs);
  EXPECT_EQ(prob_outputs.constraint_val_,
            108.0209317352670552736526587978005409240722656250000);
  EXPECT_EQ(prob_outputs.d_constraint_val_,
            -16.97102783562994332555717846844345331192016601562500);
  EXPECT_EQ(prob_outputs.d2_constraint_val_,
            23.15927346865850822155152854975312948226928710937500);
}
TEST(FakeTest, TestFake2) {
  using ::roahm::GetDeltaS;
  using ::roahm::Tick;
  const auto t0 = Tick();
  const std::string frs_in_fname{"/data/frs_bin_processed.txt"};
  std::ifstream file_in(frs_in_fname, std::ios::binary);
  const ::roahm::FrsTotal frs_in{
      ::roahm::ReadFromBinFile<::roahm::FrsTotal>(file_in)};
  file_in.close();
  const auto t1 = Tick();

  const auto& curr_frs = frs_in.megas_.at(0).au_.at(0);
  const std::array<double, 3> u0v0r0_slice_beta{0.13, 0.77, 0.04};
  const std::vector<double> mu_sigma =
      BuildMuSigma(curr_frs.cuda_info_.num_zono_);
  const std::vector<std::vector<double>> all_mu_sigmas = {mu_sigma};

  const double risk_threshold = 0.0;
  const auto t2 = Tick();
  ::roahm::risk_constraint::RiskConstraint risk_constraint{
      curr_frs, u0v0r0_slice_beta, all_mu_sigmas, 0.0};
  const auto t3 = Tick();
  auto outputs = risk_constraint.EvaluateAt(0.01, risk_threshold);
  const auto t4 = Tick();
  for (int i = 1; i < 10; ++i) {
    const auto t4_1 = Tick();
    outputs = risk_constraint.EvaluateAt(0.01, risk_threshold);
    const auto t4_2 = Tick();
    std::cout << "Time to eval cons[" + std::to_string(i) + "]: "
              << GetDeltaS(t4_2, t4_1) << "[s]" << std::endl;
  }
  const auto t5 = Tick();
  EXPECT_EQ(outputs.risk_,
            108.0209317352670552736526587978005409240722656250000);
  EXPECT_EQ(outputs.d_risk_,
            -16.97102783562994332555717846844345331192016601562500);
  EXPECT_EQ(outputs.d2_risk_,
            23.15927346865850822155152854975312948226928710937500);
  std::cout << "Time to read:         " << GetDeltaS(t1, t0) << "[s]"
            << std::endl;
  std::cout << "Time to setup cons:   " << GetDeltaS(t3, t2) << "[s]"
            << std::endl;
  std::cout << "Time to eval cons[0]: " << GetDeltaS(t4, t3) << "[s]"
            << std::endl;
  std::cout << "Time to eval cons[1]: " << GetDeltaS(t5, t4) << "[s]"
            << std::endl;
  std::cout << "Time total:           " << GetDeltaS(t5, t0) << "[s]"
            << std::endl;

  // const auto& cuda_info = frs_in.megas_.at(0).au_.at(0).cuda_info_;
  // const auto outputs = ::roahm::pre_slice::PreSlice(cuda_info,
  // u0v0r0_slice_beta, mu_sigma.data()); double p = 0.01; p = (p -
  // cuda_info.cg_p_.at(0)) / cuda_info.cg_p_.at(1); const
  // ::roahm::ProbIntegrationInputs prob_inputs  {
  //   cuda_info.num_zono_,
  //   cuda_info.grid_x0_.data(),
  //   cuda_info.grid_y0_.data(),
  //   cuda_info.grid_dx_.data(),
  //   cuda_info.grid_dy_.data(),
  //   cuda_info.block_inzono_list_.data(),
  //   cuda_info.rot_angle_.data(),
  //   mu_sigma.data(),
  //   cuda_info.cg_p_.data(),
  //   cuda_info.g_p_x_.data(),
  //   cuda_info.grid_size_,
  //   p,
  //   outputs.H1_.get(),
  //   outputs.H2_.get(),
  //   outputs.H4_.get()
  // };
  // const auto prob_outputs =
  // ::roahm::prob_integ::ProbIntegration(prob_inputs);
  // EXPECT_EQ(prob_outputs.constraint_val_,
  // 108.0209317352670552736526587978005409240722656250000);
  // EXPECT_EQ(prob_outputs.d_constraint_val_,
  // -16.97102783562994332555717846844345331192016601562500);
  // EXPECT_EQ(prob_outputs.d2_constraint_val_, 23.15927346865850822155152854975312948226928710937500);
}

TEST(FakeTest, TestFake3) {
  using ::roahm::GetDeltaS;
  using ::roahm::Tick;
  // arbitrary value
  const double risk_threshold = 91.03;

  const auto t0 = Tick();
  const std::string frs_in_fname{"/data/frs_bin_processed.txt"};
  std::ifstream file_in(frs_in_fname, std::ios::binary);
  const ::roahm::FrsTotal frs_in{
      ::roahm::ReadFromBinFile<::roahm::FrsTotal>(file_in)};
  file_in.close();
  const auto t1 = Tick();

  const auto& curr_frs = frs_in.megas_.at(0).au_.at(0);
  const std::array<double, 3> u0v0r0_slice_beta{0.13, 0.77, 0.04};
  const std::vector<double> mu_sigma =
      BuildMuSigma(curr_frs.cuda_info_.num_zono_);
  const std::vector<std::vector<double>> all_mu_sigmas = {mu_sigma};

  const auto t2 = Tick();
  ::roahm::risk_constraint::RiskConstraint risk_constraint{
      curr_frs, u0v0r0_slice_beta, all_mu_sigmas, 0.0};
  const auto t3 = Tick();
  auto outputs = risk_constraint.EvaluateAt(0.01, risk_threshold);
  const auto t4 = Tick();
  for (int i = 1; i < 10; ++i) {
    const auto t4_1 = Tick();
    outputs = risk_constraint.EvaluateAt(0.01, risk_threshold);
    const auto t4_2 = Tick();
    std::cout << "Time to eval cons[" + std::to_string(i) + "]: "
              << GetDeltaS(t4_2, t4_1) << "[s]" << std::endl;
  }
  const auto t5 = Tick();
  EXPECT_EQ(outputs.risk_,
            108.0209317352670552736526587978005409240722656250000 -
                risk_threshold);
  EXPECT_EQ(outputs.d_risk_,
            -16.97102783562994332555717846844345331192016601562500);
  EXPECT_EQ(outputs.d2_risk_,
            23.15927346865850822155152854975312948226928710937500);
  std::cout << "Time to read:         " << GetDeltaS(t1, t0) << "[s]"
            << std::endl;
  std::cout << "Time to setup cons:   " << GetDeltaS(t3, t2) << "[s]"
            << std::endl;
  std::cout << "Time to eval cons[0]: " << GetDeltaS(t4, t3) << "[s]"
            << std::endl;
  std::cout << "Time to eval cons[1]: " << GetDeltaS(t5, t4) << "[s]"
            << std::endl;
  std::cout << "Time total:           " << GetDeltaS(t5, t0) << "[s]"
            << std::endl;

  // const auto& cuda_info = frs_in.megas_.at(0).au_.at(0).cuda_info_;
  // const auto outputs = ::roahm::pre_slice::PreSlice(cuda_info,
  // u0v0r0_slice_beta, mu_sigma.data()); double p = 0.01; p = (p -
  // cuda_info.cg_p_.at(0)) / cuda_info.cg_p_.at(1); const
  // ::roahm::ProbIntegrationInputs prob_inputs  {
  //   cuda_info.num_zono_,
  //   cuda_info.grid_x0_.data(),
  //   cuda_info.grid_y0_.data(),
  //   cuda_info.grid_dx_.data(),
  //   cuda_info.grid_dy_.data(),
  //   cuda_info.block_inzono_list_.data(),
  //   cuda_info.rot_angle_.data(),
  //   mu_sigma.data(),
  //   cuda_info.cg_p_.data(),
  //   cuda_info.g_p_x_.data(),
  //   cuda_info.grid_size_,
  //   p,
  //   outputs.H1_.get(),
  //   outputs.H2_.get(),
  //   outputs.H4_.get()
  // };
  // const auto prob_outputs =
  // ::roahm::prob_integ::ProbIntegration(prob_inputs);
  // EXPECT_EQ(prob_outputs.constraint_val_,
  // 108.0209317352670552736526587978005409240722656250000);
  // EXPECT_EQ(prob_outputs.d_constraint_val_,
  // -16.97102783562994332555717846844345331192016601562500);
  // EXPECT_EQ(prob_outputs.d2_constraint_val_, 23.15927346865850822155152854975312948226928710937500);
}