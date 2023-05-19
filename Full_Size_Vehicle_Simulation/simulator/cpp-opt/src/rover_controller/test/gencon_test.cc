#include "gencon.hpp"

#include <gtest/gtest.h>

#include <chrono>

#include "dyn_obs.hpp"
#include "frs_loader.hpp"
#include "timing_util.hpp"

TEST(GenconTest, SmallTest) {
  using ::roahm::GetDeltaS;
  using ::roahm::Tick;
  const ::roahm::DynObs obs_0{1.0,
                              2.0,
                              0.1,
                              0.5,
                              5.0,
                              1.0,
                              ::roahm::DynObs::ObstacleType::kStaticBoundary};

  ::roahm::FlZonoObsSet obs_info;
  obs_info.PushObs(obs_0);
  Eigen::Matrix3Xd zono_mat;
  zono_mat.resize(3, 2);
  zono_mat << 2.1, 2.2, 4.1, 4.2, 0.1, 0.2;

  ::roahm::Vehrs vehrs;
  {
    for (int zono_idx = 0; zono_idx < 1; ++zono_idx) {
      vehrs.xy_centers_.push_back(zono_mat(0, 0));
      vehrs.xy_centers_.push_back(zono_mat(1, 0));
      vehrs.h_centers_.push_back(zono_mat(2, 0));
      for (int i = 1; i < zono_mat.cols(); ++i) {
        vehrs.zono_xy_.push_back(zono_mat(0, i));
        vehrs.zono_xy_.push_back(zono_mat(1, i));
        vehrs.zono_h_.push_back(zono_mat(2, i));
      }
      vehrs.zono_sizes_.push_back(zono_mat.cols() - 1);
      vehrs.slc_x_.push_back(0.5);
      vehrs.slc_y_.push_back(0.3);
      const double tint_0 = zono_idx * 0.01;
      const double tint_1 = tint_0 + 0.01;
      const ::roahm::Interval z_tint{tint_0, tint_1};
      vehrs.zono_time_intervals_.push_back(z_tint);
    }
  }

  const auto t0 = Tick();
  const auto cons = GenerateConstraints(vehrs, obs_info);
  const auto t1 = Tick();
  const auto cons_eig = GenerateConstraints2(vehrs, obs_info);
  const auto t2 = Tick();
  std::cout << "Time for GenerateConstraints:  " << GetDeltaS(t1, t0)
            << std::endl;
  std::cout << "Time for GenerateConstraints2: " << GetDeltaS(t2, t1)
            << std::endl;
  ASSERT_EQ(cons.num_a_elts_, cons_eig.a_mat_.size());
  ASSERT_EQ(cons.num_b_elts_, cons_eig.b_mat_.size());
  ASSERT_EQ(cons.zono_obs_sizes_.size(), cons_eig.start_idx_mat_.size());
  ASSERT_EQ(cons.zono_startpoints_.size(), cons_eig.size_mat_.size());
  ASSERT_EQ(cons_eig.start_idx_mat_.size(), cons_eig.size_mat_.size());
  std::cout << "Zono SZ" << std::endl;
  for (int i = 0; i < cons.zono_obs_sizes_.size(); ++i) {
    EXPECT_EQ(cons.zono_startpoints_[i], cons_eig.start_idx_mat_(i, 0));
    EXPECT_EQ(cons.zono_obs_sizes_[i], cons_eig.size_mat_(i, 0));
  }
  std::cout << "A Con" << std::endl;
  constexpr double kTol = 1.0e-8;
  for (int i = 0; i < cons.num_a_elts_; ++i) {
    EXPECT_NEAR(cons.a_con_arr_[i], cons_eig.a_mat_(i, 0), kTol);
  }
  std::cout << "B Con" << std::endl;
  for (int i = 0; i < cons.num_b_elts_; ++i) {
    EXPECT_NEAR(cons.b_con_arr_[i], cons_eig.b_mat_(i, 0), kTol);
  }
}

TEST(GenconTest, test1) {
  using ::roahm::GetDeltaS;
  using ::roahm::Tick;
  const ::roahm::DynObs obs_0{1.0,
                              2.0,
                              0.1,
                              0.5,
                              5.0,
                              1.0,
                              ::roahm::DynObs::ObstacleType::kStaticBoundary};
  const ::roahm::DynObs obs_1{18.0,
                              21.0,
                              -0.3,
                              0.4,
                              4.0,
                              2.0,
                              ::roahm::DynObs::ObstacleType::kStaticBoundary};

  ::roahm::FlZonoObsSet obs_info;
  for (int i = 0; i < 8; ++i) {
    obs_info.PushObs(obs_0);
    obs_info.PushObs(obs_1);
  }
  Eigen::Matrix3Xd zono_mat;
  zono_mat.resize(3, 20);
  zono_mat << 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3,
      3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8,
      4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 0.1, 0.2, 0.3,
      0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
      1.9, 2.0;

  ::roahm::Vehrs vehrs;
  {
    for (int zono_idx = 0; zono_idx < 6; ++zono_idx) {
      vehrs.xy_centers_.push_back(zono_mat(0, 0));
      vehrs.xy_centers_.push_back(zono_mat(1, 0));
      vehrs.h_centers_.push_back(zono_mat(2, 0));
      for (int i = 1; i < zono_mat.cols(); ++i) {
        vehrs.zono_xy_.push_back(zono_mat(0, i));
        vehrs.zono_xy_.push_back(zono_mat(1, i));
        vehrs.zono_h_.push_back(zono_mat(2, i));
      }
      vehrs.zono_sizes_.push_back(zono_mat.cols() - 1);
      vehrs.slc_x_.push_back(0.5);
      vehrs.slc_y_.push_back(0.3);
      const double tint_0 = zono_idx * 0.01;
      const double tint_1 = tint_0 + 0.01;
      const ::roahm::Interval z_tint{tint_0, tint_1};
      vehrs.zono_time_intervals_.push_back(z_tint);
    }
  }
  {
    const double c_x = 0.5;
    const double c_y = 0.7;
    const double c_h = 0.3;
    const double z_gx1 = 0.1;
    const double z_gx2 = -0.2;
    const double z_gx3 = -0.7;
    const double z_gy1 = 0.5;
    const double z_gy2 = 0.3;
    const double z_gy3 = -0.9;
    const double z_gh1 = 0.05;
    const double z_gh2 = -0.05;
    const double z_gh3 = -0.03;
    const ::roahm::Interval z_tint{0.0, 0.01};
    vehrs.xy_centers_.push_back(c_x);
    vehrs.xy_centers_.push_back(c_y);
    vehrs.h_centers_.push_back(c_h);
    vehrs.zono_xy_.push_back(z_gx1);
    vehrs.zono_xy_.push_back(z_gy1);
    vehrs.zono_xy_.push_back(z_gx2);
    vehrs.zono_xy_.push_back(z_gy2);
    vehrs.zono_xy_.push_back(z_gx3);
    vehrs.zono_xy_.push_back(z_gy3);
    vehrs.zono_h_.push_back(z_gh1);
    vehrs.zono_h_.push_back(z_gh2);
    vehrs.zono_h_.push_back(z_gh3);
    vehrs.zono_sizes_.push_back(3);
    vehrs.slc_x_.push_back(0.5);
    vehrs.slc_y_.push_back(0.3);
    vehrs.zono_time_intervals_.push_back(z_tint);
  }
  {
    const double c_x = 0.9;
    const double c_y = 1.1;
    const double c_h = -0.7;
    const double z_gx1 = 0.5;
    const double z_gx2 = -0.9;
    const double z_gx3 = 0.2;
    const double z_gx4 = 0.7;
    const double z_gy1 = 0.04;
    const double z_gy2 = 0.3;
    const double z_gy3 = -0.9;
    const double z_gy4 = 0.2;
    const double z_gh1 = 0.05;
    const double z_gh2 = -0.45;
    const double z_gh3 = -0.04;
    const double z_gh4 = -0.09;
    const ::roahm::Interval z_tint{0.01, 0.03};
    vehrs.xy_centers_.push_back(c_x);
    vehrs.xy_centers_.push_back(c_y);
    vehrs.h_centers_.push_back(c_h);
    vehrs.zono_xy_.push_back(z_gx1);
    vehrs.zono_xy_.push_back(z_gy1);
    vehrs.zono_xy_.push_back(z_gx2);
    vehrs.zono_xy_.push_back(z_gy2);
    vehrs.zono_xy_.push_back(z_gx3);
    vehrs.zono_xy_.push_back(z_gy3);
    vehrs.zono_xy_.push_back(z_gx4);
    vehrs.zono_xy_.push_back(z_gy4);
    vehrs.zono_h_.push_back(z_gh1);
    vehrs.zono_h_.push_back(z_gh2);
    vehrs.zono_h_.push_back(z_gh3);
    vehrs.zono_h_.push_back(z_gh4);
    vehrs.zono_sizes_.push_back(4);
    vehrs.slc_x_.push_back(0.8);
    vehrs.slc_y_.push_back(0.2);
    vehrs.zono_time_intervals_.push_back(z_tint);
  }

  const auto t0 = Tick();
  const auto cons = GenerateConstraints(vehrs, obs_info);
  const auto t1 = Tick();
  const auto cons_eig = GenerateConstraints2(vehrs, obs_info);
  const auto t2 = Tick();
  std::cout << "Time for GenerateConstraints:  " << GetDeltaS(t1, t0)
            << std::endl;
  std::cout << "Time for GenerateConstraints2: " << GetDeltaS(t2, t1)
            << std::endl;
  ASSERT_EQ(cons.num_a_elts_, cons_eig.a_mat_.size());
  ASSERT_EQ(cons.num_b_elts_, cons_eig.b_mat_.size());
  ASSERT_EQ(cons.zono_obs_sizes_.size(), cons_eig.start_idx_mat_.size());
  ASSERT_EQ(cons.zono_startpoints_.size(), cons_eig.size_mat_.size());
  ASSERT_EQ(cons_eig.start_idx_mat_.size(), cons_eig.size_mat_.size());
  std::cout << "Zono SZ" << std::endl;
  for (int i = 0; i < cons.zono_obs_sizes_.size(); ++i) {
    EXPECT_EQ(cons.zono_startpoints_[i], cons_eig.start_idx_mat_(i, 0));
    EXPECT_EQ(cons.zono_obs_sizes_[i], cons_eig.size_mat_(i, 0));
  }
  std::cout << "A Con" << std::endl;
  constexpr double kTol = 1.0e-8;
  for (int i = 0; i < cons.num_a_elts_; ++i) {
    EXPECT_NEAR(cons.a_con_arr_[i], cons_eig.a_mat_(i, 0), kTol);
  }
  std::cout << "B Con" << std::endl;
  for (int i = 0; i < cons.num_b_elts_; ++i) {
    EXPECT_NEAR(cons.b_con_arr_[i], cons_eig.b_mat_(i, 0), kTol);
  }
}