#include "triangle_gen.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <eigen3/Eigen/Dense>

#include "individual_mu_sigma.hpp"
#include "matplot/core/line_spec.h"
#include "matplot/freestanding/plot.h"
TEST(TriangleGen, TriangleGenExtraTest) {
  // Test Info:
  // * PDF is a circle with 3 sigma radius of 1, centered at (3, 2)
  // * Unsliced zonotope is a rectangle with a width of 3 and a height of 2
  // * Sliced zonotope is a rectangle with a width of 1 and a height of 1
  // * Slice generator is x-axis aligned already, with an x generator value of 2
  // This means that the PDF intersection with the unsliced should be a quarter
  // circle in the upper right of the rectangle, and the intersection with the
  // sliced should be a quarter circle in the upper right of the sliced
  // rectangle, which is a square from the origin to (1, 1).
  const double c_x_unaligned = 0.0;
  const double c_y_unaligned = 1.0;
  const std::vector<double> g_x_nonslice_unaligned{0.1, std::cos(M_PI / 4.0)};
  const std::vector<double> g_y_nonslice_unaligned{0.0, std::sin(M_PI / 4.0)};
  const double c_p = 0.0;
  const double g_p = 0.1;
  const double x_p_unaligned = 2.0;
  const double y_p_unaligned = 0.0;
  const int grid_dim = 10;
  const roahm::IndividualMuSigma mu_sigma_unaligned{3.0, 2.0, 1.0 / 9.0, 0.0,
                                                    1.0 / 9.0};
  auto ret = roahm::TriangleGenSingle(
      c_x_unaligned, c_y_unaligned, g_x_nonslice_unaligned,
      g_y_nonslice_unaligned, c_p, g_p, x_p_unaligned, y_p_unaligned, grid_dim,
      mu_sigma_unaligned, 3);
  return;
}
TEST(TriangleGen, IntersectionLineSegTest) {
  // Test Info:
  // * PDF is a circle with 3 sigma radius of 1, centered at (3, 2)
  // * Unsliced zonotope is a rectangle with a width of 3 and a height of 2
  // * Sliced zonotope is a rectangle with a width of 1 and a height of 1
  // * Slice generator is x-axis aligned already, with an x generator value of 2
  // This means that the PDF intersection with the unsliced should be a quarter
  // circle in the upper right of the rectangle, and the intersection with the
  // sliced should be a quarter circle in the upper right of the sliced
  // rectangle, which is a square from the origin to (1, 1).
  //
  // We can then compute
  const std::vector<double> g_x_nonslice_unaligned{1.0, 0.0};
  const std::vector<double> g_y_nonslice_unaligned{0.0, 1.0};
  Eigen::Matrix2Xd zono_gens;
  zono_gens.conservativeResize(2, 2);
  zono_gens << 1.0, 0.0, 0.0, 1.0;
  Eigen::Vector2d obs_center;
  obs_center << 5.0, 0.0;
  Eigen::Vector2d obs_gen;
  obs_gen << 0.0, 1.0;
  auto ret = roahm::CheckIntersectionCenteredZonoSingleGenObs(
      zono_gens, obs_center, obs_gen);
  std::cout << "Ret: " << ret << std::endl;
}

TEST(TriangleGen, TriangleGenSingleTest) {
  // Test Info:
  // * PDF is a circle with 3 sigma radius of 1, centered at (3, 2)
  // * Unsliced zonotope is a rectangle with a width of 3 and a height of 2
  // * Sliced zonotope is a rectangle with a width of 1 and a height of 1
  // * Slice generator is x-axis aligned already, with an x generator value of 2
  // This means that the PDF intersection with the unsliced should be a quarter
  // circle in the upper right of the rectangle, and the intersection with the
  // sliced should be a quarter circle in the upper right of the sliced
  // rectangle, which is a square from the origin to (1, 1).
  //
  // We can then compute
  const double c_x_unaligned = 0.0;
  const double c_y_unaligned = 0.0;
  const std::vector<double> g_x_nonslice_unaligned{1.0, 0.0};
  const std::vector<double> g_y_nonslice_unaligned{0.0, 1.0};
  const double c_p = 0.0;
  const double g_p = 0.1;
  const double x_p_unaligned = 2.0;
  const double y_p_unaligned = 2.0;
  const int grid_dim = 10;
  const roahm::IndividualMuSigma mu_sigma_unaligned{
      3.0, 1.0, 1.0 / 9.0, 1.0 * (1.0 / 9.0) / std::sqrt(2), 1.0 / 9.0};
  auto ret = roahm::TriangleGenSingle(
      c_x_unaligned, c_y_unaligned, g_x_nonslice_unaligned,
      g_y_nonslice_unaligned, c_p, g_p, x_p_unaligned, y_p_unaligned, grid_dim,
      mu_sigma_unaligned, 3);
  EXPECT_EQ(1, 1);

  // Eigen::Vector2d center{-4.0, 3.0};
  // Eigen::Matrix<double, 2, Eigen::Dynamic> generators_xy{2, 3};
  // generators_xy << 1.0, 0.0, 1.0,
  //                  0.0, 1.0, 1.0;
  // PlotZono({center, generators_xy});
  return;
}