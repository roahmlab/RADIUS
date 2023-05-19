#ifndef ROAHM_FL_ZONO_OBS_PLOT_HPP_
#define ROAHM_FL_ZONO_OBS_PLOT_HPP_
#include "dyn_obs.hpp"
#include "fl_zono_obs_set.hpp"
#include "plot_interface.hpp"
#include "point_xyh.hpp"
namespace roahm {
namespace plot_impls {
void PlotFlZonoObs(const DynObs& obs) {
  const double obs_x{obs.CenterX0()};
  const double obs_y{obs.CenterY0()};
  const double obs_h{obs.HeadingRad()};
  const double width{obs.Width()};
  const double len{obs.Length()};
  const double cos_h{std::cos(obs_h)};
  const double sin_h{std::sin(obs_h)};
  auto transform_local_pt =
      [obs_x, obs_y, cos_h,
       sin_h](const double x_local,
              const double y_local) -> std::pair<double, double> {
    const double x_rot = (x_local * cos_h + y_local * -sin_h);
    const double y_rot = (x_local * sin_h + y_local * cos_h);
    return {obs_x + x_rot, obs_y + y_rot};
  };
  const auto [x0, y0] = transform_local_pt(-len / 2.0, -width / 2.0);
  const auto [x1, y1] = transform_local_pt(+len / 2.0, -width / 2.0);
  const auto [x2, y2] = transform_local_pt(+len / 2.0, +width / 2.0);
  const auto [x3, y3] = transform_local_pt(-len / 2.0, +width / 2.0);
  std::vector<double> x{x0, x1, x2, x3, x0};
  std::vector<double> y{y0, y1, y2, y3, y0};
  plot_utils::PlotInfo plot_info{plot_utils::PlotColor{0.0, 0.0, 0.0, 0.0}, 2,
                                 false, plot_utils::PlotLineType::Continuous(),
                                 std::nullopt};
  plot_utils::PlotPolygon(x, y, plot_info);
}

void PlotFlZonoObsSet(const FlZonoObsSet& obs_set) {
  for (int i = 0; i < obs_set.GetNumObs(); ++i) {
    PlotFlZonoObs(obs_set.GetSingleDynObs(i));
    plot_utils::TurnHoldOn();
  }
}

void PlotPointXYH(const PointXYH& xyh, const double arrow_len,
                  const double tip_len, const plot_utils::PlotInfo& plot_info) {
  const double x0 = xyh.x_;
  const double y0 = xyh.y_;
  const double tip_theta = M_PI / 6.0;
  const double tip_theta_net_1 = (M_PI + xyh.h_) + tip_theta;
  const double tip_theta_net_2 = (M_PI + xyh.h_) - tip_theta;
  const double x1 = x0 + (arrow_len * std::cos(xyh.h_));
  const double y1 = y0 + (arrow_len * std::sin(xyh.h_));
  const double x2 = x1 + (tip_len * std::cos(tip_theta_net_1));
  const double y2 = y1 + (tip_len * std::sin(tip_theta_net_1));
  const double x3 = x1 + (tip_len * std::cos(tip_theta_net_2));
  const double y3 = y1 + (tip_len * std::sin(tip_theta_net_2));
  const std::vector x{x0, x1, x2, x1, x3};
  const std::vector y{y0, y1, y2, y1, y3};
  plot_utils::PlotPolygon(x, y, plot_info);
}
} // namespace plot_impls
} // namespace roahm
#endif // ROAHM_FL_ZONO_OBS_PLOT_HPP_