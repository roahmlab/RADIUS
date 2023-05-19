#include <codecvt>
#include <locale>
#include <memory>

#include "MatlabDataArray/ArrayFactory.hpp"
#include "cppmex/mexMatlabEngine.hpp"
#include "plot_interface.hpp"
#include <fmt/format.h>
#include <stdexcept>
namespace roahm {
namespace plot_utils {

std::optional<std::shared_ptr<matlab::engine::MATLABEngine>>
ConfigEngineReference(
    std::optional<std::shared_ptr<matlab::engine::MATLABEngine>> ref =
        std::nullopt) {
  static std::optional<std::shared_ptr<matlab::engine::MATLABEngine>>
      engine_ref = std::nullopt;
  if (ref.has_value() and (not engine_ref.has_value())) {
    engine_ref = ref;
  }
  return engine_ref;
}
namespace {
void ExecuteCommand(std::string cmd_in) {
  // TODO REMOVEfmt::print("[MATLAB COMMAND] {}\n", cmd_in);
  // TODO could just std::?
  std::wstring_convert<std::codecvt_utf8_utf16<char16_t>, char16_t> convert;
  std::u16string cmd = convert.from_bytes(cmd_in);
  auto engine_ptr_opt = ConfigEngineReference();
  if (engine_ptr_opt) {
    engine_ptr_opt.value()->eval(cmd);
  } else {
    throw std::runtime_error("Matlab Engine not configured yet");
  }
}
} // namespace
void PlotPolygon(std::vector<double> x, std::vector<double> y,
                 const PlotInfo& plot_info) {
  matlab::data::ArrayFactory factory;
  auto x_arr = factory.createArray<double>({x.size()});
  for (int i = 0; i < x.size(); ++i) {
    x_arr[i] = x.at(i);
  }
  auto y_arr = factory.createArray<double>({y.size()});
  for (int i = 0; i < y.size(); ++i) {
    y_arr[i] = y.at(i);
  }
  static int tmp_idx{0};
  tmp_idx++;
  const std::string x_tmp_var_name = "x_tmp______" + std::to_string(tmp_idx);
  const std::string y_tmp_var_name = "y_tmp______" + std::to_string(tmp_idx);
  auto engine_ptr_opt = ConfigEngineReference();
  if (not engine_ptr_opt) {
    throw std::runtime_error("Matlab Engine not configured yet");
  } else {
    engine_ptr_opt.value()->setVariable(x_tmp_var_name, x_arr);
    engine_ptr_opt.value()->setVariable(y_tmp_var_name, y_arr);
  }

  ExecuteCommand(fmt::format(
      "plot({}, {}, '{}', 'LineWidth', {}, 'Color', [{}, {}, {}, {}]);",
      x_tmp_var_name, y_tmp_var_name, plot_info.line_type_.GetRepresentation(),
      plot_info.line_width_, plot_info.color_.red_, plot_info.color_.green_,
      plot_info.color_.blue_, 1.0 - plot_info.color_.alpha_));
}

void OpenFigure(int fig_id) {

  ExecuteCommand(fmt::format("figure({});", fig_id));
}

void OpenSubplot(int rows, int cols, int idx) {
  ExecuteCommand(fmt::format("subplot({}, {}, {});", rows, cols, idx));
}

void TurnHoldOn() { ExecuteCommand("hold on"); }
void TurnHoldOff() { ExecuteCommand("hold off"); }
void SetTitle(std::string title) {
  ExecuteCommand(fmt::format("title(\"{}\");", title));
}
void SetAxisEqual() { ExecuteCommand("axis equal"); }
void LegendOnDisplayNamesOnly() { ExecuteCommand("legend"); }
void Show() {
  // do nothing, plots show automatically
}
} // namespace plot_utils
} // namespace roahm