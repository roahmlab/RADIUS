#include "individual_mu_sigma.hpp"

#include <gtest/gtest.h>

#include "point_xyh.hpp"

TEST(IndividualMuSigmaTest, ZeroMeanZeroBase) {
  constexpr double kEps = 1.0e-8;
  const ::roahm::IndividualMuSigma mu_sigma_global{0, 0, 1, 2, 3};
  const ::roahm::PointXYH& frame{0.0, 0.0, M_PI / 2.0};

  const auto mu_sigma_relative_no_mirror =
      mu_sigma_global.RelativeTo(frame, false);
  EXPECT_DOUBLE_EQ(mu_sigma_relative_no_mirror.mu_x_, 0.0);
  EXPECT_DOUBLE_EQ(mu_sigma_relative_no_mirror.mu_y_, 0.0);
  EXPECT_DOUBLE_EQ(mu_sigma_relative_no_mirror.sigma_1_, 3.0);
  EXPECT_DOUBLE_EQ(mu_sigma_relative_no_mirror.sigma_2_, -2.0);
  EXPECT_DOUBLE_EQ(mu_sigma_relative_no_mirror.sigma_4_, 1.0);

  const auto mu_sigma_relative_mirror = mu_sigma_global.RelativeTo(frame, true);
  EXPECT_NEAR(mu_sigma_relative_mirror.mu_x_, 0.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_mirror.mu_y_, 0.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_mirror.sigma_1_, 3.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_mirror.sigma_2_, 2.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_mirror.sigma_4_, 1.0, kEps);
}

TEST(IndividualMuSigmaTest, RotationTest) {
  constexpr double kEps = 1.0e-8;
  const ::roahm::IndividualMuSigma mu_sigma_global{0, -3, 1.0, 2.0, 4.0};
  const ::roahm::PointXYH& frame{0.0, 0.0, 3.0 * M_PI / 2.0};

  const auto mu_sigma_relative_no_mirror =
      mu_sigma_global.RelativeTo(frame, false);
  EXPECT_NEAR(mu_sigma_relative_no_mirror.mu_x_, 3.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_no_mirror.mu_y_, 0.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_no_mirror.sigma_1_, 4.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_no_mirror.sigma_2_, -2.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_no_mirror.sigma_4_, 1.0, kEps);

  const auto mu_sigma_relative_mirror = mu_sigma_global.RelativeTo(frame, true);
  EXPECT_NEAR(mu_sigma_relative_mirror.mu_x_, 3.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_mirror.mu_y_, 0.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_mirror.sigma_1_, 4.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_mirror.sigma_2_, 2.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_mirror.sigma_4_, 1.0, kEps);
}

TEST(IndividualMuSigmaTest, MultipleTest) {
  constexpr double kEps = 1.0e-8;
  const ::roahm::IndividualMuSigma mu_sigma_global{6, 7, 19, -14, 1};
  const ::roahm::PointXYH& frame{12.0, 24.0, 3.0 * M_PI / 2.0};

  const auto mu_sigma_relative_no_mirror =
      mu_sigma_global.RelativeTo(frame, false);
  EXPECT_NEAR(mu_sigma_relative_no_mirror.mu_x_, 17.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_no_mirror.mu_y_, -6.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_no_mirror.sigma_1_, 1.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_no_mirror.sigma_2_, 14.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_no_mirror.sigma_4_, 19.0, kEps);

  const auto mu_sigma_relative_mirror = mu_sigma_global.RelativeTo(frame, true);
  EXPECT_NEAR(mu_sigma_relative_mirror.mu_x_, 17.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_mirror.mu_y_, 6.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_mirror.sigma_1_, 1.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_mirror.sigma_2_, -14.0, kEps);
  EXPECT_NEAR(mu_sigma_relative_mirror.sigma_4_, 19.0, kEps);
}