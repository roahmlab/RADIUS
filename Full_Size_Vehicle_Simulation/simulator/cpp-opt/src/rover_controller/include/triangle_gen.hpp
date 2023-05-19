#ifndef ROAHM_TRIANGLE_GEN_HPP_
#define ROAHM_TRIANGLE_GEN_HPP_

#include <algorithm>
#include <array>
#include <eigen3/Eigen/Dense>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <vector>

#include "cuda_info.hpp"
#include "individual_mu_sigma.hpp"
#include "plot_interface.hpp"
#include "point_xy.hpp"
#include "point_xyh.hpp"
#include "vehrs_plot.hpp"

namespace roahm {

[[nodiscard]] constexpr static plot_utils::PlotColor SigmaColor() {
  return {0.3, 1.0, 0.0, 0.0};
}
[[nodiscard]] constexpr static plot_utils::PlotColor SigmaBoxColor() {
  return {0.0, 1.0, 0.0, 0.0};
}
[[nodiscard]] constexpr static plot_utils::PlotColor UnslicedColor() {
  return {0.5, 0.0, 0.0, 0.5};
}
[[nodiscard]] constexpr static plot_utils::PlotColor UnslicedBoxColor() {
  return {0.5, 0.0, 0.0, 0.5};
}
[[nodiscard]] constexpr static plot_utils::PlotColor SlicedColor() {
  return {0.5, 0.0, 1.0, 0.0};
}
[[nodiscard]] constexpr static plot_utils::PlotColor SlicedBoxColor() {
  return {0.0, 0.0, 1.0, 0.0};
}
[[nodiscard]] constexpr static plot_utils::PlotColor
IntersectionWithSlicedLocalFrameColor() {
  return {0.0, 1.0, 0.0, 0.5};
}

[[nodiscard]] constexpr static plot_utils::PlotInfo
BoundingBoxGeneralPlotInfo(const plot_utils::PlotColor& color,
                           std::string_view legend_entry) {
  return {color, 4, false, plot_utils::PlotLineType::Dashed(), legend_entry};
}
[[nodiscard]] constexpr static plot_utils::PlotInfo SigmaBoundingBoxPlotInfo() {
  return BoundingBoxGeneralPlotInfo(SigmaBoxColor(), "Sigma Bounding Box");
}
[[nodiscard]] constexpr static plot_utils::PlotInfo
UnslicedZonoBoundingBoxPlotInfo() {
  return BoundingBoxGeneralPlotInfo(UnslicedBoxColor(),
                                    "Unsliced Zono Bounding Box");
}

[[nodiscard]] constexpr static plot_utils::PlotInfo ObsPdfPlotInfo() {
  return {SigmaColor(), 4, true, plot_utils::PlotLineType::Continuous(),
          "Obstacle PDF"};
}
[[nodiscard]] constexpr static plot_utils::PlotInfo UnslicedZonoPlotInfo() {
  return {UnslicedColor(), 4, false, plot_utils::PlotLineType::Continuous(),
          "Unsliced Zonotope"};
}
[[nodiscard]] constexpr static plot_utils::PlotInfo SlicedZonoPlotInfo() {
  return {SlicedColor(), 4, false, plot_utils::PlotLineType::Continuous(),
          "Sliced Zonotope"};
}
[[nodiscard]] constexpr static plot_utils::PlotInfo
IntersectionWithSlicedBoundingBoxPlotInfo() {
  return {IntersectionWithSlicedLocalFrameColor(), 4, false,
          plot_utils::PlotLineType::Continuous(),
          "Intersection Sliced-Sigma (Local)"};
}

struct Triangle {
  std::array<::roahm::PointXY, 3> pts_;
};

struct BoxMinMax {
  const double x_min_;
  const double x_max_;
  const double y_min_;
  const double y_max_;
  [[nodiscard]] double Width() const {
    if (x_max_ < x_min_) {
      throw std::runtime_error("BoxMinMax x_max < x_min");
    }
    return x_max_ - x_min_;
  }
  [[nodiscard]] double Height() const {
    if (y_max_ < y_min_) {
      throw std::runtime_error("BoxMinMax y_max < y_min");
    }
    return y_max_ - y_min_;
  }
};

void PlotTriangle(const Triangle& triangle, plot_utils::PlotColor color) {
  const auto pt0 = std::get<0>(triangle.pts_);
  const auto pt1 = std::get<1>(triangle.pts_);
  const auto pt2 = std::get<2>(triangle.pts_);
  std::vector<double> x = {pt0.x_, pt1.x_, pt2.x_};
  std::vector<double> y = {pt0.y_, pt1.y_, pt2.y_};
  const bool filled = true;
  const float line_width = 2.0;
  const plot_utils::PlotLineType line_type{
      plot_utils::PlotLineType::Continuous()};
  const std::optional<std::string> legend_entry = std::nullopt;
  plot_utils::PlotInfo plot_info{color, line_width, filled, line_type,
                                 legend_entry};
  PlotPolygon(x, y, plot_info);
}

void PlotBox(const BoxMinMax& box, plot_utils::PlotInfo plot_info) {
  std::vector<double> x = {box.x_min_, box.x_max_, box.x_max_, box.x_min_,
                           box.x_min_};
  std::vector<double> y = {box.y_min_, box.y_min_, box.y_max_, box.y_max_,
                           box.y_min_};
  PlotPolygon(x, y, plot_info);
}

// TRUE iff intersection
bool CheckIntersectionCenteredZonoSingleGenObs(
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& zono_gens,
    const Eigen::Vector2d obs_center, const Eigen::Vector2d obs_gen) {
  // TODO is this c or -c?
  const auto c = obs_center;
  Eigen::Matrix2Xd big_g;
  big_g.conservativeResize(2, zono_gens.cols() + obs_gen.cols());
  big_g << zono_gens, obs_gen;
  Eigen::MatrixXd big_c;
  big_c.conservativeResizeLike(big_g);
  big_c << -big_g.row(1), big_g.row(0);
  big_c.colwise().normalize();
  big_c.transposeInPlace();
  auto delta_d =
      ((((big_c * big_g).transpose()).cwiseAbs()).colwise().sum()).transpose();
  const auto d = big_c * c;
  Eigen::MatrixXd pa;
  pa.conservativeResize(big_c.rows() * 2, big_c.cols());
  pa << big_c, -big_c;
  Eigen::MatrixXd pb;
  pb.conservativeResize(delta_d.rows() * 2, delta_d.cols());
  pb << (d + delta_d), (-d + delta_d);
  Eigen::Vector2d zero_pt;
  zero_pt.setZero();
  return 0 <= pb.minCoeff();
}

bool CheckIntersectionCenteredZono(
    const Triangle& triangle,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& generators) {
  const auto get_midpoint = [](const PointXY& p0,
                               const PointXY& p1) -> Eigen::Vector2d {
    const PointXY midpoint{(p0.x_) + ((p1.x_ - p0.x_) / 2.0),
                           (p0.y_) + ((p1.y_ - p0.y_) / 2.0)};
    Eigen::Vector2d midpoint_eig;
    midpoint_eig << midpoint.x_, midpoint.y_;
    return midpoint_eig;
  };
  const auto get_generator = [](const PointXY& p0,
                                const PointXY& p1) -> Eigen::Vector2d {
    const double dx = p1.x_ - p0.x_;
    const double dy = p1.y_ - p0.y_;
    Eigen::Vector2d generator_eig;
    generator_eig << (dx / 2.0), (dy / 2.0);
    return generator_eig;
  };
  const auto pt0 = std::get<0>(triangle.pts_);
  const auto pt1 = std::get<1>(triangle.pts_);
  const auto pt2 = std::get<2>(triangle.pts_);
  const auto midpoint01 = get_midpoint(pt0, pt1);
  const auto midpoint12 = get_midpoint(pt1, pt2);
  const auto midpoint20 = get_midpoint(pt2, pt0);
  const auto generator01 = get_generator(pt0, pt1);
  const auto generator12 = get_generator(pt1, pt2);
  const auto generator20 = get_generator(pt2, pt0);
  return CheckIntersectionCenteredZonoSingleGenObs(generators, midpoint01,
                                                   generator01) or
         CheckIntersectionCenteredZonoSingleGenObs(generators, midpoint12,
                                                   generator12) or
         CheckIntersectionCenteredZonoSingleGenObs(generators, midpoint20,
                                                   generator20);
}

void PlotMuSigma(const IndividualMuSigma& mu_sigma, const double number_stdevs,
                 std::optional<std::string> legend_entry) {
  Eigen::Matrix2d sigma_mat;
  sigma_mat << mu_sigma.sigma_1_, mu_sigma.sigma_2_, mu_sigma.sigma_2_,
      mu_sigma.sigma_4_;

  // const double pearson = mu_sigma.sigma_2_ / std::sqrt(mu_sigma.sigma_1_ *
  // mu_sigma.sigma_4_); if (pearson > 1.0 || pearson < -1.0) {
  //     throw std::runtime_error("Pearson correlation coefficient is out of
  //     range");
  // }
  // const double ellipse_radius_x = std::sqrt(1.0 + pearson);
  // const double ellipse_radius_y = std::sqrt(1.0 - pearson);
  // const double scale_x = std::sqrt(mu_sigma.sigma_1_) * number_stdevs;
  // const double scale_y = std::sqrt(mu_sigma.sigma_4_) * number_stdevs;

  Eigen::EigenSolver<Eigen::Matrix2d> eigen_solver;
  eigen_solver.compute(sigma_mat);
  std::cout << "Sigma Mat: " << sigma_mat << std::endl;
  Eigen::Vector2d eigen_values = eigen_solver.eigenvalues().real();
  Eigen::Matrix2d eigen_vectors = eigen_solver.eigenvectors().real();
  if (eigen_values.x() < 0.0 || eigen_values.y() < 0.0) {
    throw std::runtime_error("Eigen values are negative");
  }
  std::cout << "Eigen vectors [pre]: " << eigen_vectors << std::endl;
  std::cout << "Eigen values [pre]: " << eigen_values << std::endl;
  if (eigen_values.y() > eigen_values.x()) {
    eigen_values.reverseInPlace();
    eigen_vectors.rowwise().reverseInPlace();
  }
  const double major_angle =
      std::atan2(eigen_vectors(1, 0), eigen_vectors(0, 0));
  std::cout << "Eigen values: " << eigen_values << std::endl;
  std::cout << "Eigen vectors: " << eigen_vectors << std::endl;
  const double scale_x = std::sqrt(eigen_values.x()) * number_stdevs;
  const double scale_y = std::sqrt(eigen_values.y()) * number_stdevs;
  std::cout << "Scale X: " << scale_x << std::endl;
  std::cout << "Scale Y: " << scale_y << std::endl;
  std::cout << "Major angle: " << major_angle << std::endl;
  const int num_points = 100;
  Eigen::Rotation2Dd major_rotation(major_angle);
  Eigen::Matrix<double, 2, Eigen::Dynamic> points_unrotated(2, num_points);
  const double angle_step = (2.0 * M_PI) / num_points;
  for (int i = 0; i < num_points; ++i) {
    const double angle = i * angle_step;
    points_unrotated(0, i) = scale_x * std::cos(angle);
    points_unrotated(1, i) = scale_y * std::sin(angle);
  }
  const auto points_rotated =
      (major_rotation.toRotationMatrix() * points_unrotated).eval();
  Eigen::Vector2d mu_xy;
  mu_xy << mu_sigma.mu_x_, mu_sigma.mu_y_;
  std::cout << "Mu XY: \n" << mu_xy << std::endl;
  const auto points_offset = (points_rotated.colwise() + mu_xy).eval();

  std::vector<double> x;
  std::vector<double> y;
  for (int i = 0; i < points_offset.cols(); ++i) {
    x.push_back(points_offset(0, i));
    y.push_back(points_offset(1, i));
  }
  const float line_width = 4.0;
  std::cout << "plotting poly" << std::endl;
  PlotPolygon(x, y, ObsPdfPlotInfo());
  std::cout << "plotted poly" << std::endl;
}

std::optional<BoxMinMax> BoxIntersection(const BoxMinMax& box1,
                                         const BoxMinMax& box2) {
  const double x_min = std::max(box1.x_min_, box2.x_min_);
  const double x_max = std::min(box1.x_max_, box2.x_max_);
  const double y_min = std::max(box1.y_min_, box2.y_min_);
  const double y_max = std::min(box1.y_max_, box2.y_max_);
  // TODO >= should avoid measure zero intersection, but should probably set a
  // practical tolerance
  if (x_min >= x_max || y_min >= y_max) {
    return std::nullopt;
  }
  return BoxMinMax{x_min, x_max, y_min, y_max};
}

[[nodiscard]] std::optional<CudaInfo> TriangleGenSingle(
    const double c_x_unaligned, const double c_y_unaligned,
    const std::vector<double> g_x_nonslice_unaligned,
    const std::vector<double> g_y_nonslice_unaligned, const double c_p,
    const double g_p, const double g_x_p_unaligned,
    const double g_y_p_unaligned, const int grid_dim,
    const IndividualMuSigma mu_sigma_unaligned, const double num_stdevs) {
  // Plot Initial Data

  // auto fig = matplot::figure();
  // fig->size(800, 1500);
  constexpr int kNumSubplotRows = 3;
  constexpr int kNumSubplotCols = 1;
  constexpr int kUnalignedSubplotIdx = 0;
  constexpr int kAlignedSubplotIdx = 1;
  constexpr int kSlicedSubplotIdx = 2;
  {
    constexpr auto is_valid_subplot = [](const int idx) {
      return idx >= 0 and idx < (kNumSubplotRows * kNumSubplotCols);
    };
    static_assert(is_valid_subplot(kUnalignedSubplotIdx));
    static_assert(is_valid_subplot(kAlignedSubplotIdx));
    static_assert(is_valid_subplot(kSlicedSubplotIdx));
  }

  plot_utils::OpenSubplot(kNumSubplotRows, kNumSubplotCols,
                          kUnalignedSubplotIdx);
  {
    auto g_x_unsliced_unaligned = g_x_nonslice_unaligned;
    g_x_unsliced_unaligned.push_back(g_x_p_unaligned);
    auto g_y_unsliced_unaligned = g_y_nonslice_unaligned;
    g_y_unsliced_unaligned.push_back(g_y_p_unaligned);
    PlotZono({c_x_unaligned, c_y_unaligned, g_x_unsliced_unaligned,
              g_y_unsliced_unaligned},
             UnslicedZonoPlotInfo());
    plot_utils::TurnHoldOn();
    PlotZono({c_x_unaligned, c_y_unaligned, g_x_nonslice_unaligned,
              g_y_nonslice_unaligned},
             SlicedZonoPlotInfo());
    PlotMuSigma(mu_sigma_unaligned, num_stdevs, "Mu Sigma");
  }
  plot_utils::TurnHoldOff();
  plot_utils::SetTitle("Unaligned");
  plot_utils::SetAxisEqual();
  plot_utils::LegendOnDisplayNamesOnly();

  if (g_x_nonslice_unaligned.size() != g_y_nonslice_unaligned.size()) {
    throw std::runtime_error(
        "g_x_nonslice_unaligned.size() != g_y_nonslice_unaligned.size()");
  }
  const double slice_gen_angle = std::atan2(g_y_p_unaligned, g_x_p_unaligned);
  const Eigen::Rotation2Dd rot_mat_align(-slice_gen_angle);
  Eigen::Matrix<double, 2, Eigen::Dynamic> g_nonslice_unaligned;
  g_nonslice_unaligned.resize(2, g_x_nonslice_unaligned.size());
  for (int col = 0; col < g_x_nonslice_unaligned.size(); ++col) {
    g_nonslice_unaligned.col(col) << g_x_nonslice_unaligned.at(col),
        g_y_nonslice_unaligned.at(col);
  }
  const auto g_nonslice_aligned =
      (rot_mat_align.toRotationMatrix() * g_nonslice_unaligned).eval();
  std::cout << "g_nonslice_unaligned: " << g_nonslice_unaligned << std::endl;
  std::cout << "slice_gen_angle: " << slice_gen_angle << std::endl;
  std::cout << "g_nonslice_aligned: " << g_nonslice_aligned << std::endl;

  double x_abs{0.0};
  double y_abs{0.0};
  for (int col = 0; col < g_nonslice_aligned.cols(); ++col) {
    x_abs += std::abs(g_nonslice_aligned(0, col));
    y_abs += std::abs(g_nonslice_aligned(1, col));
  }

  const auto g_xy_p_aligned =
      (rot_mat_align * Eigen::Vector2d(g_x_p_unaligned, g_y_p_unaligned))
          .eval();
  const double g_x_p_aligned = g_xy_p_aligned(0);
  const double g_y_p_aligned = g_xy_p_aligned(1);
  std::cout << "g_x_p_aligned: " << g_x_p_aligned << std::endl;
  std::cout << "g_y_p_aligned: " << g_y_p_aligned << std::endl;
  const double sliced_zono_box_width = 2 * x_abs;
  const double sliced_zono_box_height = 2 * y_abs;
  const double unsliced_zono_box_width =
      sliced_zono_box_width + 2 * std::abs(g_x_p_aligned);
  const double unsliced_zono_box_height =
      sliced_zono_box_height + 2 * std::abs(g_y_p_aligned);
  const auto mu_sigma_aligned = mu_sigma_unaligned.RelativeTo(
      PointXYH{c_x_unaligned, c_y_unaligned, slice_gen_angle}, false);
  std::cout << "Sliced Zono Box Width: " << sliced_zono_box_width << std::endl;
  std::cout << "Sliced Zono Box Height : " << sliced_zono_box_height
            << std::endl;
  std::cout << "Unsliced Zono Box Width: " << unsliced_zono_box_width
            << std::endl;
  std::cout << "Unsliced Zono Box Height : " << unsliced_zono_box_height
            << std::endl;
  Eigen::Matrix2d sigma_aligned;
  sigma_aligned << mu_sigma_aligned.sigma_1_, mu_sigma_aligned.sigma_2_,
      mu_sigma_aligned.sigma_2_, mu_sigma_aligned.sigma_4_;
  Eigen::EigenSolver<Eigen::Matrix2d> eigen_solver;
  eigen_solver.compute(sigma_aligned);
  // Variances along the major and minor axes
  Eigen::Vector2d eigen_values = eigen_solver.eigenvalues().real();
  if (eigen_values.x() < 0.0 || eigen_values.y() < 0.0) {
    throw std::runtime_error("Eigen values are negative");
  }
  // The standard deviations along the major and minor axes
  Eigen::Vector2d eigen_values_sqrt = eigen_values.cwiseSqrt();
  Eigen::Matrix2d eigen_vectors = eigen_solver.eigenvectors().real();
  std::cout << "Eigen Values: " << eigen_values << std::endl;
  std::cout << "Eigen Vectors: " << eigen_vectors << std::endl;
  constexpr double kSigmaToUse = 3;
  // TODO check this is correct
  const double sigma_box_width =
      kSigmaToUse * (std::abs(eigen_values_sqrt(0) * eigen_vectors(0, 0)) +
                     std::abs(eigen_values_sqrt(1) * eigen_vectors(0, 1)));
  const double sigma_box_height =
      kSigmaToUse * (std::abs(eigen_values_sqrt(0) * eigen_vectors(1, 0)) +
                     std::abs(eigen_values_sqrt(1) * eigen_vectors(1, 1)));
  const double sigma_box_x_min_aligned =
      mu_sigma_aligned.mu_x_ - sigma_box_width;
  const double sigma_box_y_min_aligned =
      mu_sigma_aligned.mu_y_ - sigma_box_height;
  const double sigma_box_x_max_aligned =
      mu_sigma_aligned.mu_x_ + sigma_box_width;
  const double sigma_box_y_max_aligned =
      mu_sigma_aligned.mu_y_ + sigma_box_height;
  std::cout << "Sigma Box Width: " << sigma_box_width << std::endl;
  std::cout << "Sigma Box Height: " << sigma_box_height << std::endl;
  std::cout << "Sigma Box X Min: " << sigma_box_x_min_aligned << std::endl;
  std::cout << "Sigma Box Y Min: " << sigma_box_y_min_aligned << std::endl;
  std::cout << "Sigma Box X Max: " << sigma_box_x_max_aligned << std::endl;
  std::cout << "Sigma Box Y Max: " << sigma_box_y_max_aligned << std::endl;

  // Find intersection of unsliced zono box with sigma box
  const BoxMinMax sigma_box_aligned{
      sigma_box_x_min_aligned, sigma_box_x_max_aligned, sigma_box_y_min_aligned,
      sigma_box_y_max_aligned};
  const BoxMinMax unsliced_zono_box_aligned{
      -unsliced_zono_box_width / 2, unsliced_zono_box_width / 2,
      -unsliced_zono_box_height / 2, unsliced_zono_box_height / 2};
  const BoxMinMax sliced_zono_box_aligned{
      -sliced_zono_box_width / 2, sliced_zono_box_width / 2,
      -sliced_zono_box_height / 2, sliced_zono_box_height / 2};
  const auto intersection_unsliced_with_sigma =
      BoxIntersection(sigma_box_aligned, unsliced_zono_box_aligned);
  std::cout << "unsliced_zono_box_aligned.x: "
            << unsliced_zono_box_aligned.x_min_ << ", "
            << unsliced_zono_box_aligned.x_max_ << std::endl;
  std::cout << "unsliced_zono_box_aligned.y: "
            << unsliced_zono_box_aligned.y_min_ << ", "
            << unsliced_zono_box_aligned.y_max_ << std::endl;

  // Plot Aligned Data
  plot_utils::OpenSubplot(kNumSubplotRows, kNumSubplotCols, kAlignedSubplotIdx);
  {
    Eigen::Matrix<double, 2, Eigen::Dynamic> g_aligned_with_slice_gen;
    g_aligned_with_slice_gen.resize(2, g_nonslice_aligned.cols() + 1);
    g_aligned_with_slice_gen << g_nonslice_aligned, g_xy_p_aligned;
    PlotZono({{0.0, 0.0}, g_aligned_with_slice_gen}, UnslicedZonoPlotInfo());
    plot_utils::TurnHoldOn();
    PlotZono({{0.0, 0.0}, g_nonslice_aligned}, SlicedZonoPlotInfo());
    PlotMuSigma(mu_sigma_aligned, num_stdevs, "Mu Sigma");
    PlotBox(sigma_box_aligned, SigmaBoundingBoxPlotInfo());
    PlotBox(unsliced_zono_box_aligned, UnslicedZonoBoundingBoxPlotInfo());
  }
  plot_utils::TurnHoldOff();
  plot_utils::SetTitle("Aligned");
  plot_utils::SetAxisEqual();
  plot_utils::LegendOnDisplayNamesOnly();

  if (not intersection_unsliced_with_sigma) {
    std::cout << "No intersection between sigma box and unsliced zono box"
              << std::endl;
    plot_utils::Show();
    return std::nullopt;
  }
  const double slice_gen_length = std::abs(g_x_p_aligned);
  const double sigma_box_x_aligned_buffered_width =
      std::abs(sigma_box_width) + std::abs(slice_gen_length);
  const double sigma_box_x_min_aligned_buffered =
      mu_sigma_aligned.mu_x_ - sigma_box_x_aligned_buffered_width;
  const double sigma_box_x_max_aligned_buffered =
      mu_sigma_aligned.mu_x_ + sigma_box_x_aligned_buffered_width;
  const BoxMinMax sigma_box_aligned_buffered{
      sigma_box_x_min_aligned_buffered, sigma_box_x_max_aligned_buffered,
      sigma_box_y_min_aligned, sigma_box_y_max_aligned};
  const auto intersection_sliced_with_sigma_opt =
      BoxIntersection(sigma_box_aligned_buffered, sliced_zono_box_aligned);
  if (not intersection_sliced_with_sigma_opt) {
    plot_utils::Show();
    throw std::runtime_error(
        "No intersection between sigma box and sliced zono box but "
        "intersection with unsliced zono box!");
    return std::nullopt;
  }

  const auto intersection_sliced_with_sigma =
      intersection_sliced_with_sigma_opt.value();
  std::cout << "Intersection Sliced With Sigma X Min: "
            << intersection_sliced_with_sigma.x_min_ << std::endl;
  std::cout << "Intersection Sliced With Sigma X Max: "
            << intersection_sliced_with_sigma.x_max_ << std::endl;
  std::cout << "Intersection Sliced With Sigma Y Min: "
            << intersection_sliced_with_sigma.y_min_ << std::endl;
  std::cout << "Intersection Sliced With Sigma Y Max: "
            << intersection_sliced_with_sigma.y_max_ << std::endl;

  std::vector<Triangle> triangles;
  const double grid_dx =
      intersection_sliced_with_sigma.Width() / static_cast<double>(grid_dim);
  const double grid_dy =
      intersection_sliced_with_sigma.Height() / static_cast<double>(grid_dim);
  const double grid_x0_local = intersection_sliced_with_sigma.x_min_;
  const double grid_y0_local = intersection_sliced_with_sigma.y_min_;
  // Starting from LL
  for (int row = 0; row < grid_dim; ++row) {
    for (int col = 0; col < grid_dim; ++col) {
      const double r_x = ((col + 1) * grid_dx) + grid_x0_local;
      const double l_x = (col * grid_dx) + grid_x0_local;
      const double l_y = (row * grid_dy) + grid_y0_local;
      const double u_y = ((row + 1) * grid_dy) + grid_y0_local;

      const ::roahm::PointXY ll_pt{l_x, l_y};
      const ::roahm::PointXY ul_pt{l_x, u_y};
      const ::roahm::PointXY lr_pt{r_x, l_y};
      const ::roahm::PointXY ur_pt{r_x, u_y};

      // TODO should these be LL/UR or LR/UL??
      const auto triangle_ll =
          Triangle{std::array<::roahm::PointXY, 3>{ll_pt, ul_pt, lr_pt}};
      const auto triangle_ur =
          Triangle{std::array<::roahm::PointXY, 3>{ul_pt, ur_pt, lr_pt}};
      triangles.push_back(triangle_ll);
      triangles.push_back(triangle_ur);
    }
  }

  // Plot Aligned Data
  plot_utils::OpenFigure(2);
  // matplot::subplot(kNumSubplotRows, kNumSubplotCols, kSlicedSubplotIdx);
  {
    Eigen::Matrix<double, 2, Eigen::Dynamic> g_aligned_with_slice_gen;
    g_aligned_with_slice_gen.resize(2, g_nonslice_aligned.cols() + 1);
    g_aligned_with_slice_gen << g_nonslice_aligned, g_xy_p_aligned;
    PlotZono({{0.0, 0.0}, g_aligned_with_slice_gen}, UnslicedZonoPlotInfo());
    plot_utils::TurnHoldOn();
    PlotZono({{0.0, 0.0}, g_nonslice_aligned}, SlicedZonoPlotInfo());
    PlotMuSigma(mu_sigma_aligned, num_stdevs, "Mu Sigma");
    PlotBox(sigma_box_aligned, SigmaBoundingBoxPlotInfo());
    PlotBox(sigma_box_aligned_buffered, SigmaBoundingBoxPlotInfo());
    PlotBox(intersection_sliced_with_sigma,
            IntersectionWithSlicedBoundingBoxPlotInfo());
    for (int i = 0; i < triangles.size(); i += 2) {
      const auto triangle_ll = triangles.at(i + 0);
      const auto triangle_ur = triangles.at(i + 1);
      const bool ll_intersects =
          CheckIntersectionCenteredZono(triangle_ll, g_nonslice_aligned);
      const bool ur_intersects =
          CheckIntersectionCenteredZono(triangle_ur, g_nonslice_aligned);
      const plot_utils::PlotColor intersect_color =
          plot_utils::PlotColor{0.9, 1.0, 0.0, 0.0};
      const plot_utils::PlotColor no_intersect_color =
          plot_utils::PlotColor{0.9, 0.0, 1.0, 0.0};
      const plot_utils::PlotColor ll_color =
          ll_intersects ? intersect_color : no_intersect_color;
      const plot_utils::PlotColor ur_color =
          ur_intersects ? intersect_color : no_intersect_color;
      PlotTriangle(triangle_ll, ll_color);
      PlotTriangle(triangle_ur, ur_color);
    }
  }
  plot_utils::TurnHoldOff();
  plot_utils::SetTitle("Sliced Grid");
  plot_utils::SetAxisEqual();
  // matplot::xlim({-10, 10});
  // matplot::ylim({-10, 10});
  plot_utils::LegendOnDisplayNamesOnly();
  plot_utils::Show();
  // TODO actually generate a return value

  return std::nullopt;
}

} // namespace roahm
#endif // ROAHM_TRIANGLE_GEN_HPP_