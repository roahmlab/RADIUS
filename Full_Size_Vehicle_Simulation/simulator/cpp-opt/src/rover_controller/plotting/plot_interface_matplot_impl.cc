#include "plot_interface.hpp"

#include <matplot/matplot.h>

#include "matplot/core/legend.h"
#include "matplot/util/common.h"
namespace roahm {
namespace plot_utils {

void PlotPolygon(std::vector<double> x, std::vector<double> y,
                 const PlotInfo& plot_info) {
  const auto& filled = plot_info.filled_;
  const auto& color = plot_info.color_;
  const auto& legend_entry = plot_info.legend_entry_;
  const auto& line_width = plot_info.line_width_;
  const std::string_view line_style = plot_info.line_type_.GetRepresentation();
  auto lh = matplot::polygon(x, y, "g");
  lh->fill(filled).color({color.alpha_, color.red_, color.green_, color.blue_});
  if (legend_entry.has_value()) {
    lh->display_name(legend_entry.value());
  }
  lh->line_style(line_style);
  lh->line_width(line_width);
}

void OpenFigure(int fig_id) {
  auto fig = matplot::figure(fig_id);
  fig->ioff();
  fig->quiet_mode(true);
}

void OpenSubplot(int rows, int cols, int idx) {
  matplot::subplot(rows, cols, idx);
}

void TurnHoldOn() { matplot::hold(matplot::on); }
void TurnHoldOff() { matplot::hold(matplot::off); }
void SetTitle(std::string title) { matplot::title(title); }

void SetAxisEqual() { matplot::axis(matplot::equal); }

void LegendOnDisplayNamesOnly() {
  matplot::legend(std::vector<std::string>{""});
}

void Show() { matplot::show(); }

} // namespace plot_utils

} // namespace roahm