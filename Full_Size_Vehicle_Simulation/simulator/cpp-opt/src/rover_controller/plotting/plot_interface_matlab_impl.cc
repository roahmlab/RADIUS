#include "MatlabDataArray/CharArray.hpp"
#include "mex.hpp"
#include "plot_interface.hpp"
#include <memory>
#include <optional>
#include <string>
#include <tuple>

class MatlabPlotter {
  std::optional<std::shared_ptr<matlab::engine::MATLABEngine>> engine_ptr_opt_;
  std::optional<matlab::data::Array> axes_handle_opt_;
  matlab::data::ArrayFactory factory_;
  std::optional<matlab::data::Array> GetAxHandle() {
    if (not axes_handle_opt_ and engine_ptr_opt_) {
      Figure();
    }

    if (axes_handle_opt_) {
      return axes_handle_opt_.value();
    }
    return std::nullopt;
  }

public:
  MatlabPlotter(std::shared_ptr<matlab::engine::MATLABEngine> engine_ptr)
      : engine_ptr_opt_{engine_ptr},
        axes_handle_opt_{std::nullopt}, factory_{} {}

  void Figure(int fig_no = 99) {
    if (not engine_ptr_opt_) {
      return;
    }
    matlab::data::TypedArray<double> fig_no_arg =
        factory_.createScalar<double>(fig_no);
    matlab::data::Array fig_ax_out;
    std::vector<matlab::data::Array> args{fig_no_arg};
    axes_handle_opt_ = engine_ptr_opt_.value()->feval(u"blue_fig", args);
  }
  void ClearAxes() {
    if (not engine_ptr_opt_) {
      return;
    }
    if (axes_handle_opt_.has_value()) {
      engine_ptr_opt_.value()->feval(u"cla", 0, {axes_handle_opt_.value()});
    }
  }

  void HoldOn() {
    if (not engine_ptr_opt_) {
      return;
    }
    const auto hold_on = factory_.createCharArray("on");
    if (axes_handle_opt_.has_value()) {
      engine_ptr_opt_.value()->feval(u"hold", 0,
                                     {axes_handle_opt_.value(), hold_on});
    }
  }
  void HoldOff() {
    if (not engine_ptr_opt_) {
      return;
    }
    const auto hold_off = factory_.createCharArray("off");
    if (axes_handle_opt_) {
      engine_ptr_opt_.value()->feval(u"hold", 0,
                                     {axes_handle_opt_.value(), hold_off});
    }
  }

  void Plot(const std::vector<double>& x, const std::vector<double>& y,
            const roahm::plot_utils::PlotInfo& plot_info) {
    const int n = x.size();
    if (not engine_ptr_opt_) {
      return;
    }
    matlab::data::TypedArray<double> x_vals =
        factory_.createArray<double>({static_cast<unsigned long>(n)});
    matlab::data::TypedArray<double> y_vals =
        factory_.createArray<double>({static_cast<unsigned long>(n)});
    for (int i = 0; i < n; ++i) {
      x_vals[i] = x.at(i);
      y_vals[i] = y.at(i);
    }
    const matlab::data::CharArray line_width_key =
        factory_.createCharArray("LineWidth");
    const matlab::data::TypedArray<double> line_width_val =
        factory_.createScalar<double>(plot_info.line_width_);
    const matlab::data::CharArray line_type_char_arr = factory_.createCharArray(
        std::string(plot_info.line_type_.GetRepresentation()));
    if (engine_ptr_opt_ and axes_handle_opt_) {
      // TODO actually implement filled
      if (true or not plot_info.filled_) {
        const matlab::data::CharArray color_key =
            factory_.createCharArray("Color");
        matlab::data::TypedArray<double> color_rgba_arr =
            factory_.createArray<double>({4});
        color_rgba_arr[0] = plot_info.color_.red_;
        color_rgba_arr[1] = plot_info.color_.green_;
        color_rgba_arr[2] = plot_info.color_.blue_;
        color_rgba_arr[3] = 1.0 - plot_info.color_.alpha_;
        engine_ptr_opt_.value()->feval(
            u"plot",
            {axes_handle_opt_.value(), x_vals, y_vals, line_type_char_arr,
             line_width_key, line_width_val, color_key, color_rgba_arr});
      } else {
        const matlab::data::CharArray line_style_key =
            factory_.createCharArray("LineStyle");
        const matlab::data::CharArray face_alpha_key =
            factory_.createCharArray("FaceAlpha");
        matlab::data::TypedArray<double> color_rgb_arr =
            factory_.createArray<double>({3});
        color_rgb_arr[0] = plot_info.color_.red_;
        color_rgb_arr[1] = plot_info.color_.green_;
        color_rgb_arr[2] = plot_info.color_.blue_;
        matlab::data::TypedArray<double> alpha_val =
            factory_.createScalar<double>(plot_info.line_width_);
        alpha_val[0] = 1.0 - plot_info.color_.alpha_;
        // TODO
        // FIXME
        engine_ptr_opt_.value()->feval(
            u"patch", 0,
            {axes_handle_opt_.value(), x_vals, y_vals, color_rgb_arr,
             line_width_key, line_width_val, line_style_key, line_type_char_arr,
             face_alpha_key, alpha_val});
      }
    }
  }
  void ConfigEngine(std::shared_ptr<matlab::engine::MATLABEngine> engine_ptr) {
    engine_ptr_opt_ = engine_ptr;
  }

  void SetAxisEqual() {
    if (engine_ptr_opt_ and axes_handle_opt_) {
      const auto equal_str = factory_.createCharArray("equal");
      engine_ptr_opt_.value()->feval(u"axis", 0,
                                     {axes_handle_opt_.value(), equal_str});
    }
  }
};

namespace roahm {
namespace plot_utils {

std::optional<std::shared_ptr<MatlabPlotter>>
GetPlotter(std::optional<std::shared_ptr<matlab::engine::MATLABEngine>>
               engine_ptr_opt = std::nullopt) {
  static std::optional<std::shared_ptr<MatlabPlotter>> plotter{};
  if (engine_ptr_opt) {
    plotter = std::make_shared<MatlabPlotter>(engine_ptr_opt.value());
  }
  return plotter;
}
void SetupEngine(std::shared_ptr<matlab::engine::MATLABEngine> engine_ptr) {
  std::ignore = GetPlotter(engine_ptr);
}

void PlotPolygon(std::vector<double> x, std::vector<double> y,
                 const PlotInfo& plot_info) {
  if (auto plotter_opt = GetPlotter(); plotter_opt.has_value()) {
    plotter_opt.value()->Plot(x, y, plot_info);
  }
}

void OpenFigure(int fig_id) {
  if (auto plotter_opt = GetPlotter(); plotter_opt.has_value()) {
    plotter_opt.value()->Figure(fig_id);
  }
}

void OpenSubplot(int rows, int cols, int idx) {
  // TODO
  return;
}

void TurnHoldOn() {
  if (auto plotter_opt = GetPlotter(); plotter_opt.has_value()) {
    plotter_opt.value()->HoldOn();
  }
}
void TurnHoldOff() {
  if (auto plotter_opt = GetPlotter(); plotter_opt.has_value()) {
    plotter_opt.value()->HoldOff();
  }
}
void SetTitle(std::string title) {
  // TODO
  // if (auto plotter_opt = GetPlotter(); plotter_opt.has_value()) {
  //  plotter_opt.value()->SetTitle(title);
  //}
}
void SetAxisEqual() {
  if (auto plotter_opt = GetPlotter(); plotter_opt.has_value()) {
    plotter_opt.value()->SetAxisEqual();
  }
}
void ClearAxes() {
  if (auto plotter_opt = GetPlotter(); plotter_opt.has_value()) {
    plotter_opt.value()->ClearAxes();
  }
}
void LegendOnDisplayNamesOnly() {
  // TODO
}
void Show() { return; }
} // namespace plot_utils
} // namespace roahm