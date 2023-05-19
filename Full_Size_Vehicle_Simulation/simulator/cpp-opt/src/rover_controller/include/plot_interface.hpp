#ifndef ROAHM_PLOT_UTILS_HPP_
#define ROAHM_PLOT_UTILS_HPP_
#include <optional>
#include <string>
#include <vector>

namespace roahm {
namespace plot_utils {

struct PlotColor {
  float alpha_;
  float red_;
  float green_;
  float blue_;
};

struct PlotLineType {
private:
  const std::string_view representation_;
  PlotLineType() = delete;
  constexpr explicit PlotLineType(std::string_view representation)
      : representation_{representation} {}

public:
  [[nodiscard]] constexpr static PlotLineType Dashed() {
    return PlotLineType{"--"};
  }
  [[nodiscard]] constexpr static PlotLineType Continuous() {
    return PlotLineType{"-"};
  };
  [[nodiscard]] constexpr std::string_view GetRepresentation() const {
    return representation_;
  }
};

struct PlotInfo {
  plot_utils::PlotColor color_;
  float line_width_;
  bool filled_;
  PlotLineType line_type_;
  std::optional<std::string_view> legend_entry_;
};

void PlotPolygon(std::vector<double> x, std::vector<double> y,
                 const PlotInfo& plot_info);

void OpenFigure(int fig_id);

void OpenSubplot(int rows, int cols, int idx);

void TurnHoldOn();
void TurnHoldOff();
void SetTitle(std::string title);
void SetAxisEqual();
void LegendOnDisplayNamesOnly();
void Show();
} // namespace plot_utils
} // namespace roahm
#endif // ROAHM_PLOT_UTILS_HPP_