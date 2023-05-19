#include "dyn_obs.hpp"

#include <gtest/gtest.h>

bool InPlusMinusPi(double x) { return (x >= -M_PI) and (x <= M_PI); }

bool AnglesNearEqual(const double theta1_rad, const double theta2_rad) {
  const double kEpsilon = 1e-6;
  const double cos_1 = std::cos(theta1_rad);
  const double sin_1 = std::sin(theta1_rad);
  const double cos_2 = std::cos(theta2_rad);
  const double sin_2 = std::sin(theta2_rad);
  return (std::abs(cos_1 - cos_2) < kEpsilon) &&
         (std::abs(sin_1 - sin_2) < kEpsilon);
}

TEST(DynObs, Construction) {
  const double x = 1.0;
  const double y = 2.0;
  const double v = 5.0;
  const double l = 6.0;
  const double w = 7.0;
  for (double h = -18.0; h < 18.0; h += 3.0) {
    const ::roahm::DynObs obs{
        x, y, h, v, l, w, ::roahm::DynObs::ObstacleType::kStaticBoundary};
    EXPECT_EQ(obs.CenterX0(), x);
    EXPECT_EQ(obs.CenterY0(), y);
    EXPECT_TRUE(AnglesNearEqual(obs.HeadingRad(), h));
    EXPECT_TRUE(InPlusMinusPi(obs.HeadingRad()));
    EXPECT_EQ(obs.Velocity(), v);
    EXPECT_EQ(obs.Length(), l);
    EXPECT_EQ(obs.Width(), w);
  }
}
TEST(DynObs, UnmirroredHeadingOnly) {
  const double x = std::sqrt(2) * 2.0;
  const double y = std::sqrt(2) * 2.0;
  const double h = 4.0;
  const double v = 5.0;
  const double l = 6.0;
  const double w = 7.0;
  constexpr double kTol = 1.0e-10;

  {
    const ::roahm::DynObs obs_base{
        x, y, h, v, l, w, ::roahm::DynObs::ObstacleType::kStaticBoundary};
    const ::roahm::PointXYH origin{0.0, 0.0, M_PI / 4.0};
    const auto obs_rel = obs_base.RelativeTo(origin, false);
    EXPECT_NEAR(obs_rel.CenterX0(), 4.0, kTol);
    EXPECT_NEAR(obs_rel.CenterY0(), 0.0, kTol);
    EXPECT_TRUE(AnglesNearEqual(obs_rel.HeadingRad(), h - origin.h_));
    EXPECT_TRUE(InPlusMinusPi(obs_rel.HeadingRad()));
    EXPECT_EQ(obs_rel.Velocity(), v);
    EXPECT_EQ(obs_rel.Length(), l);
    EXPECT_EQ(obs_rel.Width(), w);
  }
  {
    const ::roahm::DynObs obs_base{
        x, y, h, v, l, w, ::roahm::DynObs::ObstacleType::kStaticBoundary};
    const ::roahm::PointXYH origin{0.0, 0.0, -M_PI / 4.0};
    const auto obs_rel = obs_base.RelativeTo(origin, false);
    EXPECT_NEAR(obs_rel.CenterX0(), 0.0, kTol);
    EXPECT_NEAR(obs_rel.CenterY0(), 4.0, kTol);
    EXPECT_TRUE(AnglesNearEqual(obs_rel.HeadingRad(), h - origin.h_));
    EXPECT_TRUE(InPlusMinusPi(obs_rel.HeadingRad()));
    EXPECT_EQ(obs_rel.Velocity(), v);
    EXPECT_EQ(obs_rel.Length(), l);
    EXPECT_EQ(obs_rel.Width(), w);
  }
  {
    const ::roahm::DynObs obs_base{
        x, y, h, v, l, w, ::roahm::DynObs::ObstacleType::kStaticBoundary};
    const ::roahm::PointXYH origin{0.0, 0.0, 5.0 * M_PI / 4.0};
    const auto obs_rel = obs_base.RelativeTo(origin, false);
    EXPECT_NEAR(obs_rel.CenterX0(), -4.0, kTol);
    EXPECT_NEAR(obs_rel.CenterY0(), 0.0, kTol);
    EXPECT_TRUE(AnglesNearEqual(obs_rel.HeadingRad(), h - origin.h_));
    EXPECT_TRUE(InPlusMinusPi(obs_rel.HeadingRad()));
    EXPECT_EQ(obs_rel.Velocity(), v);
    EXPECT_EQ(obs_rel.Length(), l);
    EXPECT_EQ(obs_rel.Width(), w);
  }
}

TEST(DynObs, MirroredHeadingOnly) {
  const double x = std::sqrt(2) * 2.0;
  const double y = std::sqrt(2) * 2.0;
  const double h = 4.0;
  const double v = 5.0;
  const double l = 6.0;
  const double w = 7.0;
  constexpr double kTol = 1.0e-10;

  {
    const ::roahm::DynObs obs_base{
        x, y, h, v, l, w, ::roahm::DynObs::ObstacleType::kStaticBoundary};
    const ::roahm::PointXYH origin{0.0, 0.0, M_PI / 4.0};
    const auto obs_rel = obs_base.RelativeTo(origin, true);
    EXPECT_NEAR(obs_rel.CenterX0(), 4.0, kTol);
    EXPECT_NEAR(obs_rel.CenterY0(), 0.0, kTol);
    EXPECT_TRUE(AnglesNearEqual(obs_rel.HeadingRad(), -(h - origin.h_)));
    EXPECT_TRUE(InPlusMinusPi(obs_rel.HeadingRad()));
    EXPECT_EQ(obs_rel.Velocity(), v);
    EXPECT_EQ(obs_rel.Length(), l);
    EXPECT_EQ(obs_rel.Width(), w);
  }
  {
    const ::roahm::DynObs obs_base{
        x, y, h, v, l, w, ::roahm::DynObs::ObstacleType::kStaticBoundary};
    const ::roahm::PointXYH origin{0.0, 0.0, -M_PI / 4.0};
    const auto obs_rel = obs_base.RelativeTo(origin, true);
    EXPECT_NEAR(obs_rel.CenterX0(), 0.0, kTol);
    EXPECT_NEAR(obs_rel.CenterY0(), -4.0, kTol);
    EXPECT_TRUE(AnglesNearEqual(obs_rel.HeadingRad(), -(h - origin.h_)));
    EXPECT_TRUE(InPlusMinusPi(obs_rel.HeadingRad()));
    EXPECT_EQ(obs_rel.Velocity(), v);
    EXPECT_EQ(obs_rel.Length(), l);
    EXPECT_EQ(obs_rel.Width(), w);
  }
  {
    const ::roahm::DynObs obs_base{
        x, y, h, v, l, w, ::roahm::DynObs::ObstacleType::kStaticBoundary};
    const ::roahm::PointXYH origin{0.0, 0.0, 5.0 * M_PI / 4.0};
    const auto obs_rel = obs_base.RelativeTo(origin, true);
    EXPECT_NEAR(obs_rel.CenterX0(), -4.0, kTol);
    EXPECT_NEAR(obs_rel.CenterY0(), 0.0, kTol);
    EXPECT_TRUE(AnglesNearEqual(obs_rel.HeadingRad(), -(h - origin.h_)));
    EXPECT_TRUE(InPlusMinusPi(obs_rel.HeadingRad()));
    EXPECT_EQ(obs_rel.Velocity(), v);
    EXPECT_EQ(obs_rel.Length(), l);
    EXPECT_EQ(obs_rel.Width(), w);
  }
}

TEST(DynObs, BothHeadingAndTranslation) {
  const double x = 3.0;
  const double y = 5.0;
  const double h = 9.0;
  const double v = 5.0;
  const double l = 6.0;
  const double w = 7.0;
  constexpr double kTol = 1.0e-10;

  {
    const ::roahm::DynObs obs_base{
        x, y, h, v, l, w, ::roahm::DynObs::ObstacleType::kStaticBoundary};
    const ::roahm::PointXYH origin{1.0, 2.0, M_PI / 2.0};
    const auto obs_rel = obs_base.RelativeTo(origin, false);
    EXPECT_NEAR(obs_rel.CenterX0(), 3.0, kTol);
    EXPECT_NEAR(obs_rel.CenterY0(), -2.0, kTol);
    EXPECT_TRUE(AnglesNearEqual(obs_rel.HeadingRad(), h - origin.h_));
    EXPECT_TRUE(InPlusMinusPi(obs_rel.HeadingRad()));
    EXPECT_EQ(obs_rel.Velocity(), v);
    EXPECT_EQ(obs_rel.Length(), l);
    EXPECT_EQ(obs_rel.Width(), w);
  }
  {
    const ::roahm::DynObs obs_base{
        x, y, h, v, l, w, ::roahm::DynObs::ObstacleType::kStaticBoundary};
    const ::roahm::PointXYH origin{1.0, 2.0, M_PI / 2.0};
    const auto obs_rel = obs_base.RelativeTo(origin, true);
    EXPECT_NEAR(obs_rel.CenterX0(), 3.0, kTol);
    EXPECT_NEAR(obs_rel.CenterY0(), 2.0, kTol);
    EXPECT_TRUE(AnglesNearEqual(obs_rel.HeadingRad(), -(h - origin.h_)));
    EXPECT_TRUE(InPlusMinusPi(obs_rel.HeadingRad()));
    EXPECT_EQ(obs_rel.Velocity(), v);
    EXPECT_EQ(obs_rel.Length(), l);
    EXPECT_EQ(obs_rel.Width(), w);
  }
}