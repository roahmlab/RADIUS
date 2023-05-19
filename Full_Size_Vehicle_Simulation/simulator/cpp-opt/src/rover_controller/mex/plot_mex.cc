#include "MatlabDataArray/ArrayFactory.hpp"
#include "MatlabDataArray/ArrayType.hpp"
#include "MatlabDataArray/MDArray.hpp"
#include "MatlabDataArray/String.hpp"
#include "MatlabDataArray/TypedArray.hpp"
#include "mexAdapter.hpp"
#include "plot_interface.hpp"

/*
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

  void HoldOn() {
    if (not engine_ptr_opt_) {
      return;
    }
    const auto hold_on = factory_.createCharArray("on");
    if (axes_handle_opt_) {
      engine_ptr_opt_.value()->feval(u"hold", 0,
                                     {axes_handle_opt_.value(), hold_on});
    }
  }

  void Plot(const std::vector<double>& x, const std::vector<double>& y) {
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
    if (engine_ptr_opt_ and axes_handle_opt_) {
      engine_ptr_opt_.value()->feval(
          u"plot", {axes_handle_opt_.value(), x_vals, y_vals});
    }
  }
};
*/

namespace roahm::plot_utils {
void SetupEngine(std::shared_ptr<matlab::engine::MATLABEngine>);
}

class MexFunction : public matlab::mex::Function {
private:
public:
  MexFunction() {}
  void operator()(matlab::mex::ArgumentList outputs,
                  matlab::mex::ArgumentList inputs) {
    std::shared_ptr<matlab::engine::MATLABEngine> engine_ptr = getEngine();
    matlab::data::ArrayFactory factory{};
    matlab::data::TypedArray<double> fig_no = factory.createScalar<double>(13);
    matlab::data::Array fig_ax_out;
    std::vector<matlab::data::Array> args{fig_no};

    // matlab::data::TypedArray<matlab::data::String> hold_on_args{hold_on_arg};
    // MatlabPlotter plotter{engine_ptr};
    ::roahm::plot_utils::SetupEngine(engine_ptr);
    ::roahm::plot_utils::OpenFigure(25);
    std::vector<double> x_vals;
    std::vector<double> y_vals;
    for (int i = 0; i < 3; ++i) {
      x_vals.push_back(i);
      y_vals.push_back(2 * i * 3 * i);
    }
    roahm::plot_utils::PlotInfo plot_inf{
        roahm::plot_utils::PlotColor{},
        4,
        false,
        roahm::plot_utils::PlotLineType::Continuous(),
        {}};
    ::roahm::plot_utils::PlotPolygon(x_vals, y_vals, plot_inf);
    ::roahm::plot_utils::TurnHoldOn();
    x_vals.clear();
    y_vals.clear();
    for (int i = 0; i < 3; ++i) {
      x_vals.push_back(i);
      y_vals.push_back(2 + i);
    }
    ::roahm::plot_utils::PlotPolygon(x_vals, y_vals, plot_inf);
    ::roahm::plot_utils::SetAxisEqual();
    /*
    fig_ax_out = engine_ptr->feval(u"blue_fig", args);
    const int N = 3;
    matlab::data::TypedArray<double> x_vals = factory.createArray<double>({N});
    matlab::data::TypedArray<double> y_vals = factory.createArray<double>({N});
    for (int i = 0; i < N; ++i) {
      x_vals[i] = 2 * i;
      y_vals[i] = (2 * i) * (3 * i);
    }
    engine_ptr->feval(u"plot", {fig_ax_out, x_vals, y_vals});
    for (int i = 0; i < N; ++i) { x_vals[i] = i;
      y_vals[i] = i;
    }
    engine_ptr->feval(u"plot", {fig_ax_out, x_vals, y_vals});
    */
  }
};