
#include "waypoint_cost.hpp"

#include <gtest/gtest.h>

#include "point_xyh.hpp"
#include "waypoint.hpp"

/*
TEST(WaypointCostTest, ZeroDesZeroCenterNoEpsilon) {
  const ::roahm::WaypointLocalMirrorTakenIntoAccount des_wp{
      ::roahm::PointXYH{0.0, 0.0, 0.0}};
  const double sliceable_center = 1.0;
  const double sliceable_generator_val = 2.0;
  const double x_center = 0.0;
  const double y_center = 0.0;
  const double h_center = 0.0;

  // Perfect square on x^2+y^2 for very simple check
  const double sliceable_x_generator = 3.0;
  const double sliceable_y_generator = 4.0;

  const double sliceable_h_generator = 9.0;
  const double xy_weight = 3.0;
  const double h_weight = 4.0;
  // Note that due to epsilon = 0, some evaluation points
  // could risk division by zero
  const double xy_sqrt_epsilon = 0.0;
  const double h_sqrt_epsilon = 0.0;

  const ::roahm::WaypointCost wp_cost(
      des_wp, sliceable_center, sliceable_generator_val, x_center, y_center,
      h_center, sliceable_x_generator, sliceable_y_generator,
      sliceable_h_generator, xy_weight, h_weight, xy_sqrt_epsilon,
      h_sqrt_epsilon);
  // Could change these to EXPECT_NEAR
  {
    const auto ret = wp_cost.EvaluateAt(3.0);
    EXPECT_DOUBLE_EQ(ret.cost_, 51.0);
    EXPECT_DOUBLE_EQ(ret.cost_gradient_wrt_param_, 25.5);
    EXPECT_DOUBLE_EQ(ret.cost_hessian_wrt_param_, 0.0);
  }
  {
    const auto ret = wp_cost.EvaluateAt(-3.0);
    EXPECT_DOUBLE_EQ(ret.cost_, 102.0);
    EXPECT_DOUBLE_EQ(ret.cost_gradient_wrt_param_, -25.5);
    EXPECT_DOUBLE_EQ(ret.cost_hessian_wrt_param_, 0.0);
  }
  {
    const auto ret = wp_cost.EvaluateAt(0.0);
    EXPECT_DOUBLE_EQ(ret.cost_, 25.5);
    EXPECT_DOUBLE_EQ(ret.cost_gradient_wrt_param_, -25.5);
    EXPECT_DOUBLE_EQ(ret.cost_hessian_wrt_param_, 0.0);
  }
  {
    const auto ret = wp_cost.EvaluateAt(1.5);
    EXPECT_DOUBLE_EQ(ret.cost_, 12.75);
    EXPECT_DOUBLE_EQ(ret.cost_gradient_wrt_param_, 25.5);
    EXPECT_DOUBLE_EQ(ret.cost_hessian_wrt_param_, 0.0);
  }
}
*/

TEST(WaypointCostTest, NumericalDiffTest) {
  const roahm::WaypointLocalMirrorTakenIntoAccount des_wp{
      ::roahm::PointXYH{19.7, -14.3, 0.8}};
  const double sliceable_center = 1.2;
  const double sliceable_generator_val = -0.3;
  const double x_center = 9.3;
  const double y_center = -18.4;
  const double h_center = 0.2;
  const double sliceable_x_generator = 3.9;
  const double sliceable_y_generator = -5.7;

  const double sliceable_h_generator = 1.1;
  const double x_weight = 3.7;
  const double y_weight = 19.0;
  const double h_weight = 4.8;
  // Note that due to epsilon = 0, some evaluation points
  // could risk division by zero
  const double x_sqrt_epsilon = 3.0 * 1.0e-2;
  const double y_sqrt_epsilon = 4.0 * 1.0e-2;
  const double h_sqrt_epsilon = 5.0 * 1.0e-3;

  const ::roahm::WaypointCost wp_cost{des_wp,
                                      sliceable_center,
                                      sliceable_generator_val,
                                      x_center,
                                      y_center,
                                      h_center,
                                      sliceable_x_generator,
                                      sliceable_y_generator,
                                      sliceable_h_generator,
                                      x_weight,
                                      y_weight,
                                      h_weight,
                                      x_sqrt_epsilon,
                                      y_sqrt_epsilon,
                                      h_sqrt_epsilon};

  const double kEps{1.0e-10};
  const auto exact_out = wp_cost.EvaluateAt(1.5);
  const auto p0 = wp_cost.EvaluateAt(1.5 + kEps);
  const double numerical_grad{(p0.cost_ - exact_out.cost_) / kEps};
  const double numerical_grad2{
      (p0.cost_gradient_wrt_param_ - exact_out.cost_gradient_wrt_param_) /
      (kEps)};
  std::cout << "Numerical Grad: " << numerical_grad << std::endl;
  std::cout << "Exact Grad: " << exact_out.cost_gradient_wrt_param_
            << std::endl;
  std::cout << "Numerical Grad 2: " << numerical_grad2 << std::endl;
  std::cout << "Exact Grad 2: " << exact_out.cost_hessian_wrt_param_
            << std::endl;
}

/*
TEST(WaypointCostTest, FullTest) {
  const roahm::WaypointLocalMirrorTakenIntoAccount des_wp{
      ::roahm::PointXYH{19.7, -14.3, 0.8}};
  const double sliceable_center = 1.2;
  const double sliceable_generator_val = -0.3;
  const double x_center = 9.3;
  const double y_center = -18.4;
  const double h_center = 0.2;
  const double sliceable_x_generator = 3.9;
  const double sliceable_y_generator = -5.7;

  const double sliceable_h_generator = 1.1;
  const double xy_weight = 3.7;
  const double h_weight = 4.8;
  // Note that due to epsilon = 0, some evaluation points
  // could risk division by zero
  const double xy_sqrt_epsilon = 3.0 * 1.0e-2;
  const double h_sqrt_epsilon = 5.0 * 1.0e-3;

  const ::roahm::WaypointCost wp_cost(
      des_wp, sliceable_center, sliceable_generator_val, x_center, y_center,
      h_center, sliceable_x_generator, sliceable_y_generator,
      sliceable_h_generator, xy_weight, h_weight, xy_sqrt_epsilon,
      h_sqrt_epsilon);
  // Could change these to EXPECT_NEAR
  constexpr double kEpsTol = 1.0e-7;
  {
    const auto ret = wp_cost.EvaluateAt(1.5);
    EXPECT_NEAR(ret.cost_, 61.4110713992, kEpsTol);
    EXPECT_NEAR(ret.cost_gradient_wrt_param_, 73.1994392985, kEpsTol);
    EXPECT_NEAR(ret.cost_hessian_wrt_param_, 78.2472748095, kEpsTol);
  }
  {
    const auto ret = wp_cost.EvaluateAt(1.0);
    EXPECT_NEAR(ret.cost_, 41.8060958211, kEpsTol);
    EXPECT_NEAR(ret.cost_gradient_wrt_param_, -31.7774828005, kEpsTol);
    EXPECT_NEAR(ret.cost_hessian_wrt_param_, 264.066987807, kEpsTol);
  }
  {
    const auto ret = wp_cost.EvaluateAt(-1.0);
    EXPECT_NEAR(ret.cost_, 218.536202539, kEpsTol);
    EXPECT_NEAR(ret.cost_gradient_wrt_param_, -100.678313221, kEpsTol);
    EXPECT_NEAR(ret.cost_hessian_wrt_param_, 1.93604161082, kEpsTol);
  }
}
*/
