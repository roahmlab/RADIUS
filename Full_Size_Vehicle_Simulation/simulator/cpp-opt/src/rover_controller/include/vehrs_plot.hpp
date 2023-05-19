#ifndef ROAHM_VEHRS_PLOT_HPP_
#define ROAHM_VEHRS_PLOT_HPP_
#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <numeric>

#include "plot_interface.hpp"
#include "vehrs.hpp"
namespace roahm {

struct EigenZono {
  Eigen::Vector2d center_;
  Eigen::Matrix2Xd generators_xy_;
  EigenZono(const Eigen::Vector2d& center,
            const Eigen::Matrix<double, 2, Eigen::Dynamic>& generators_xy)
      : center_{center}, generators_xy_{generators_xy} {}

  EigenZono(const double c_x, const double c_y, const std::vector<double>& g_x,
            const std::vector<double>& g_y) {
    center_ = Eigen::Vector2d(c_x, c_y);
    generators_xy_ = Eigen::Matrix<double, 2, Eigen::Dynamic>(2, g_x.size());
    for (int i = 0; i < g_x.size(); ++i) {
      generators_xy_(0, i) = g_x.at(i);
      generators_xy_(1, i) = g_y.at(i);
    }
  }
};

void GetZonoPoints(const EigenZono& zono, std::vector<double>& x_vals,
                   std::vector<double>& y_vals) {
  auto center = zono.center_;
  auto generators_xy = zono.generators_xy_;
  const auto xy_max = generators_xy.cwiseAbs().rowwise().sum();

  // Flip any generators that are in the negative y direction
  std::transform(generators_xy.colwise().begin(), generators_xy.colwise().end(),
                 generators_xy.colwise().begin(),
                 [](const Eigen::Vector2d& v) -> Eigen::Vector2d {
                   if (v.y() < 0) {
                     return -v;
                   }
                   return v;
                 });

  const auto handle_angle = [](const double theta_rad) -> double {
    if (theta_rad < 0) {
      return theta_rad + 2.0 * M_PI;
    }
    return theta_rad;
  };

  const auto get_angle = [handle_angle](const Eigen::Vector2d& v) -> double {
    return handle_angle(atan2(v.y(), v.x()));
  };
  std::vector<std::size_t> col_idxs_sorted_by_angle(generators_xy.cols());
  std::iota(col_idxs_sorted_by_angle.begin(), col_idxs_sorted_by_angle.end(),
            0);
  std::sort(col_idxs_sorted_by_angle.begin(), col_idxs_sorted_by_angle.end(),
            [&generators_xy, get_angle](const std::size_t idx0,
                                        std::size_t idx1) -> bool {
              return get_angle(generators_xy.col(idx0)) <
                     get_angle(generators_xy.col(idx1));
            });
  Eigen::Matrix<double, 2, Eigen::Dynamic> sorted_cumsum_gens{
      2, generators_xy.cols() + 1};
  sorted_cumsum_gens.setZero();
  for (int i = 1; i < sorted_cumsum_gens.cols(); ++i) {
    const int idx_to_pull = col_idxs_sorted_by_angle[i - 1];
    const Eigen::Vector2d& gen_to_pull = 2 * generators_xy.col(idx_to_pull);
    const Eigen::Vector2d& prev_cumsum = sorted_cumsum_gens.col(i - 1);
    sorted_cumsum_gens.col(i) = gen_to_pull + prev_cumsum;
  }
  const double max_p_x = sorted_cumsum_gens.row(0).maxCoeff();
  const double cumsum_offset_x = xy_max.x() - max_p_x;
  const double cumsum_offset_y = -xy_max.y();
  Eigen::Vector2d cumsum_offset;
  cumsum_offset << cumsum_offset_x, cumsum_offset_y;
  sorted_cumsum_gens.colwise() += cumsum_offset;
  const Eigen::Vector2d cumsum_mirror_offset =
      sorted_cumsum_gens.col(0) +
      sorted_cumsum_gens.col(sorted_cumsum_gens.cols() - 1);
  Eigen::Matrix<double, 2, Eigen::Dynamic> output_points{
      2, sorted_cumsum_gens.cols() * 2};
  output_points.block(0, 0, 2, sorted_cumsum_gens.cols()) = sorted_cumsum_gens;
  const auto sub_block = -sorted_cumsum_gens;
  const auto mirror_chunk = (sub_block.colwise() + cumsum_mirror_offset).eval();
  auto insert_block = output_points.block(0, sorted_cumsum_gens.cols(), 2,
                                          sorted_cumsum_gens.cols());
  insert_block = mirror_chunk;
  output_points = output_points.colwise() + center;
  for (int i = 0; i < output_points.cols(); ++i) {
    x_vals.push_back(output_points(0, i));
    y_vals.push_back(output_points(1, i));
  }
}

void PlotZono(const EigenZono& zono, plot_utils::PlotInfo plot_info) {
  std::vector<double> x_vals;
  std::vector<double> y_vals;
  GetZonoPoints(zono, x_vals, y_vals);
  PlotPolygon(x_vals, y_vals, plot_info);
}

void PlotVehrs(const Vehrs& vehrs, const PointXYH& base_location,
               const bool mirror,
               std::optional<plot_utils::PlotInfo> plot_info_optional) {
  long gens_offset{0};
  std::vector<double> x_vals;
  std::vector<double> y_vals;
  x_vals.reserve(10000);
  y_vals.reserve(10000);

  for (int i = 0; i < vehrs.GetNumZonos(); ++i) {
    const long num_gens = static_cast<long>(vehrs.zono_sizes_.at(i));
    if (i % 1 == 0) {
      Eigen::Matrix2Xd gens;
      gens.conservativeResize(2, num_gens);
      for (int gen_idx = 0; gen_idx < num_gens; ++gen_idx) {
        gens(0, gen_idx) = vehrs.zono_xy_.at(2 * (gen_idx + gens_offset));
        gens(1, gen_idx) = vehrs.zono_xy_.at(2 * (gen_idx + gens_offset) + 1);
      }
      Eigen::Vector2d center;
      center << vehrs.xy_centers_.at(i * 2), vehrs.xy_centers_.at(i * 2 + 1);
      EigenZono zono{center, gens};
      GetZonoPoints(zono, x_vals, y_vals);
      x_vals.push_back(std::numeric_limits<double>::quiet_NaN());
      y_vals.push_back(std::numeric_limits<double>::quiet_NaN());
      // PlotZono(zono, plot_info);
    }
    plot_utils::TurnHoldOn();
    gens_offset += num_gens;
  }
  {
    const double cos_h = std::cos(base_location.h_);
    const double sin_h = std::sin(base_location.h_);
    const double x_offset = base_location.x_;
    const double y_offset = base_location.y_;
    const int mirror_mult = mirror ? -1.0 : 1.0;
    for (int i = 0; i < x_vals.size(); ++i) {
      if (mirror) {
      }
      const double x_in = x_vals.at(i);
      const double y_in = mirror_mult * y_vals.at(i);
      const double x_out = (x_in * cos_h + y_in * -sin_h) + x_offset;
      const double y_out = (x_in * sin_h + y_in * cos_h) + y_offset;
      x_vals.at(i) = x_out;
      y_vals.at(i) = y_out;
    }
  }
  const plot_utils::PlotInfo plot_info_default{
      plot_utils::PlotColor{0.01, 0.0, 0.0, 1.0}, 1, false,
      plot_utils::PlotLineType::Continuous()};
  plot_utils::PlotPolygon(x_vals, y_vals,
                          plot_info_optional.value_or(plot_info_default));
}

} // namespace roahm
#endif