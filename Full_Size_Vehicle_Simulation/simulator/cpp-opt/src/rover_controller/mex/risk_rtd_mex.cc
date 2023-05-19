#include <cstdio>
#include <exception>
#include <fmt/format.h>

#include <array>
#include <chrono>
#include <memory>
#include <omp.h>
#include <optional>
#include <stdexcept>
#include <string>

#include "MatlabDataArray/ArrayFactory.hpp"
#include "MatlabDataArray/TypedArray.hpp"
#include "MatlabDataArray/detail/ArrayFactoryHelpers.hpp"
#include "comparison_methods.hpp"
#include "dyn_obs.hpp"
#include "fl_zono_obs_plot.hpp"
#include "fl_zono_obs_set.hpp"
#include "frs_loader.hpp"
#include "frs_select_info.hpp"
#include "frs_total.hpp"
#include "gencon.hpp"
#include "ipopt_string_utils.hpp"
#include "manu_type.hpp"
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "mu_sigma_multi.hpp"
#include "plot_interface.hpp"
#include "point_xyh.hpp"
#include "risk_rtd.hpp"
#include "risk_rtd_ipopt_problem.hpp"
#include "ros/time.h"
#include "rover_state.hpp"
#include "timing_util.hpp"
#include "vehrs_plot.hpp"
#include "waypoint.hpp"

namespace roahm::plot_utils {
void SetupEngine(std::shared_ptr<matlab::engine::MATLABEngine> engine_ptr);
void ClearAxes();
// TODO REMOVE
// std::optional<std::shared_ptr<matlab::engine::MATLABEngine>>
// ConfigEngineReference(
//     std::optional<std::shared_ptr<matlab::engine::MATLABEngine>> ref =
//         std::nullopt);
} // namespace roahm::plot_utils

int64_t LoadIntFromArray(const matlab::data::Array& arr) {
  if (arr.getNumberOfElements() != 1) {
    throw std::runtime_error(
        fmt::format("LoadIntFromArray expected 1 element but got {}",
                    arr.getNumberOfElements()));
  }
  const auto arr_t = arr.getType();
  if ((arr_t != matlab::data::ArrayType::INT64) and
      (arr_t != matlab::data::ArrayType::INT32) and
      (arr_t != matlab::data::ArrayType::INT16) and
      (arr_t != matlab::data::ArrayType::INT8)) {
    throw std::runtime_error("LoadIntFromArray Expected integer type");
  }
  // TODO does this work on multidim?
  return arr[0];
}
std::vector<double> LoadVecDblFromArray(const matlab::data::Array& arr) {
  const auto arr_t = arr.getType();
  if (arr_t != matlab::data::ArrayType::DOUBLE) {
    throw std::runtime_error("LoadVecDblFromArray Expected type double");
  }
  const matlab::data::TypedArray<double>& arr_dbl{arr};
  std::vector<double> ret;
  ret.reserve(arr_dbl.getNumberOfElements());
  for (const double& val : arr_dbl) {
    ret.push_back(val);
  }
  return ret;
}

template <typename T>
::roahm::WaypointGlobalNoMirror
LoadWaypointFromMatlab(const T& waypoint_mat) noexcept(false) {
  constexpr std::size_t kExpectedNumElts = 3;
  const auto num_wp_elts = waypoint_mat.getNumberOfElements();
  if (num_wp_elts != kExpectedNumElts) {
    throw std::runtime_error(
        "Error: waypoint has " + std::to_string(num_wp_elts) +
        " components, instead of " + std::to_string(kExpectedNumElts));
  }

  const auto dims = waypoint_mat.getDimensions();
  if (dims.empty()) {
    throw std::runtime_error("Error: waypoint has empty dimension list!");
  }
  if (dims.at(0) != kExpectedNumElts) {
    throw std::runtime_error(
        "Error: waypoint has " + std::to_string(dims.at(0)) +
        " dimensions, instead of " + std::to_string(kExpectedNumElts));
  }
  const double x = waypoint_mat[0];
  const double y = waypoint_mat[1];
  const double h = waypoint_mat[2];
  return roahm::WaypointGlobalNoMirror{::roahm::PointXYH{x, y, h}};
}

::roahm::RoverState
LoadRoverStateFromMatlab(const matlab::data::Array& state_var) noexcept(false) {
  if (state_var.getType() != matlab::data::ArrayType::DOUBLE) {
    throw std::runtime_error("Error: rover state type is not DOUBLE!\n");
  }
  constexpr std::size_t kExpectedNumElts = 6;
  const auto num_state_elts = state_var.getNumberOfElements();
  if (num_state_elts != kExpectedNumElts) {
    throw std::runtime_error(
        "Error: state has " + std::to_string(num_state_elts) +
        " components, instead of " + std::to_string(kExpectedNumElts));
  }

  const auto dims = state_var.getDimensions();
  if (dims.empty()) {
    throw std::runtime_error("Error: state has empty dimension list!");
  }
  if (dims.at(0) != kExpectedNumElts) {
    throw std::runtime_error("Error: state has " + std::to_string(dims.at(0)) +
                             " dimensions, instead of " +
                             std::to_string(kExpectedNumElts));
  }
  const double x = state_var[0];
  const double y = state_var[1];
  const double h = state_var[2];
  const double u = state_var[3];
  const double v = state_var[4];
  const double r = state_var[5];
  const double w = u;
  return ::roahm::RoverState{x, y, h, u, v, r, w};
}

template <typename T>
std::vector<double>
ValidateAndLoadMuSigmaVec(const T& mu_sigma_mat) noexcept(false) {
  if (mu_sigma_mat.getType() != matlab::data::ArrayType::DOUBLE) {
    throw std::runtime_error("Error: mu sigma type is not DOUBLE!\n");
  }
  const int mu_sigma_num_dims = mu_sigma_mat.getDimensions().size();
  std::cout << "Mu Sigma num dims: " << mu_sigma_num_dims << std::endl;

  if (mu_sigma_num_dims != 2) {
    throw std::runtime_error("Error: mu sigma is not 2D, it has " +
                             std::to_string(mu_sigma_num_dims) +
                             " dimensions!");
  }

  const int mu_sigma_num_rows = mu_sigma_mat.getDimensions()[0];
  if (mu_sigma_num_rows != 1) {
    throw std::runtime_error("Error: mu sigma has " +
                             std::to_string(mu_sigma_num_rows) +
                             " rows, expected 1!");
  }
  constexpr int kMinMuSigmaCols = 5;
  constexpr int kMuSigmaColMultiples = 5;
  const int num_mu_sigma_cols = mu_sigma_mat.getDimensions()[1];
  if (num_mu_sigma_cols < 5) {
    throw std::runtime_error(
        "Error: mu sigma has " + std::to_string(num_mu_sigma_cols) +
        " cols, expected at least " + std::to_string(kMinMuSigmaCols) + "!");
  }
  if (num_mu_sigma_cols % kMuSigmaColMultiples != 0) {
    throw std::runtime_error("Error: mu sigma has " +
                             std::to_string(num_mu_sigma_cols) +
                             " cols, expected a multiple of " +
                             std::to_string(kMuSigmaColMultiples) + "!");
  }

  std::vector<double> ret;
  ret.reserve(num_mu_sigma_cols);
  for (int i = 0; i < num_mu_sigma_cols; ++i) {
    ret.push_back(mu_sigma_mat[i]);
  }
  return ret;
}

std::vector<roahm::MuSigmaMulti>
ProcessAlwaysRiskyInput(const matlab::data::Array& arr) {
  if (arr.getType() != matlab::data::ArrayType::STRUCT) {
    throw std::runtime_error("Always risky obstacles do not have type STRUCT");
  } else {
    // TODO REMOVE
    fmt::print("[DBG] arr has type struct\n");
  }

  constexpr int kNumExpectedDims = 2;
  const int num_dims = arr.getDimensions().size();
  fmt::print("Number of dimension: {}\n", num_dims);
  if (num_dims != kNumExpectedDims) {
    throw std::runtime_error(fmt::format("Expected {} dimensions but got {}",
                                         kNumExpectedDims, num_dims));
  }

  const int num_mu_sigma_sets = arr.getDimensions()[0];
  fmt::print("Dim[0]: {}\n", num_mu_sigma_sets);
  if (num_mu_sigma_sets <= 0) {
    return {};
  }

  constexpr int kNumExpectedMuSigmaDurations = 16;
  const int num_mu_sigma_durations = arr.getDimensions()[1];
  fmt::print("Dim[1]: {}\n", num_mu_sigma_durations);
  if (num_mu_sigma_durations != kNumExpectedMuSigmaDurations) {
    throw std::runtime_error(fmt::format("Expected {} cols but got {}",
                                         kNumExpectedMuSigmaDurations,
                                         num_mu_sigma_durations));
  }

  const matlab::data::StructArray& mat_struct_arr{arr};
  const auto& field_names_mat_str = mat_struct_arr.getFieldNames();
  constexpr std::string_view kMuSigmaDurationMsFieldName{
      "mu_sigma_duration_ms"};
  constexpr std::string_view kMuSigmaDtMsFieldName{"mu_sigma_dt_ms"};
  constexpr std::string_view kMuSigmaDataFieldName{"mu_sigma_data"};
  constexpr int kNumExpectedFields = 3;
  constexpr std::array<std::string_view, kNumExpectedFields>
      kExpectedFieldNames{kMuSigmaDtMsFieldName, kMuSigmaDurationMsFieldName,
                          kMuSigmaDataFieldName};
  if (kNumExpectedFields != mat_struct_arr.getNumberOfFields()) {
    throw std::runtime_error(fmt::format("Expected {} fields but got {}",
                                         kNumExpectedFields,
                                         mat_struct_arr.getNumberOfFields()));
  }
  int found_fields = 0;
  for (const auto& field_name_mat_str : mat_struct_arr.getFieldNames()) {
    for (const auto& field_name_expceted : kExpectedFieldNames) {
      if (field_name_expceted == std::string(field_name_mat_str)) {
        ++found_fields;
        break;
      }
    }
  }
  if (found_fields != kNumExpectedFields) {
    throw std::runtime_error(fmt::format("Found {} expected fields out of {}",
                                         found_fields, kNumExpectedFields));
  }

  std::vector<roahm::MuSigmaMulti> mu_sigma_multis;
  mu_sigma_multis.reserve(num_mu_sigma_sets);
  for (int mu_sigma_set_idx = 0; mu_sigma_set_idx < num_mu_sigma_sets;
       ++mu_sigma_set_idx) {
    std::array<std::vector<double>, kNumExpectedMuSigmaDurations> mu_sigma_data;
    std::optional<int64_t> mu_sigma_dt_ms{std::nullopt};
    std::optional<int64_t> mu_sigma_duration_ms{std::nullopt};
    for (int mu_sigma_duration_idx = 0;
         mu_sigma_duration_idx < kNumExpectedMuSigmaDurations;
         ++mu_sigma_duration_idx) {
      const matlab::data::Struct& single_mu_sigma =
          mat_struct_arr[mu_sigma_set_idx][mu_sigma_duration_idx];
      const int64_t single_mu_sigma_dt_ms{LoadIntFromArray(
          single_mu_sigma[std::string{kMuSigmaDtMsFieldName}])};
      if (mu_sigma_dt_ms.has_value()) {
        if (single_mu_sigma_dt_ms != mu_sigma_dt_ms) {
          throw std::runtime_error(
              fmt::format("Previously loaded dt = {}ms, curent is {}ms",
                          mu_sigma_dt_ms.value(), single_mu_sigma_dt_ms));
        }
      } else {
        mu_sigma_dt_ms = single_mu_sigma_dt_ms;
      }
      const int64_t single_mu_sigma_duration_ms{LoadIntFromArray(
          single_mu_sigma[std::string{kMuSigmaDurationMsFieldName}])};
      // TODO could do check here, but lets not for now
      mu_sigma_duration_ms = single_mu_sigma_duration_ms;
      mu_sigma_data.at(mu_sigma_duration_idx) = LoadVecDblFromArray(
          single_mu_sigma[std::string{kMuSigmaDataFieldName}]);
    }
    roahm::TimeMicroseconds dt_us{
        roahm::TimeMicroseconds::FromMilliseconds(mu_sigma_dt_ms.value())};
    mu_sigma_multis.emplace_back(mu_sigma_data, dt_us);
  }

  fmt::print("Return Size: {}\n", mu_sigma_multis.size());
  return mu_sigma_multis;
}

::roahm::FlZonoObsSet ProcessDynObs(const matlab::data::Array& data,
                                    ::roahm::DynObs::ObstacleType obs_type) {
  if (data.isEmpty() or data.getNumberOfElements() == 0) {
    return {};
  }
  constexpr int kNumExpectedDims = 2;
  if (data.getDimensions().size() != kNumExpectedDims) {
    throw std::runtime_error(fmt::format("Expected {} dimensions but had {}",
                                         kNumExpectedDims,
                                         data.getDimensions().size()));
  }

  constexpr int kNumExpectedRows = 6;
  const auto num_rows = data.getDimensions()[0];
  if (num_rows != kNumExpectedRows) {
    throw std::runtime_error(
        fmt::format("Expected {} rows but got {}", kNumExpectedRows, num_rows));
  }

  const auto num_cols = data.getDimensions()[1];
  ::roahm::FlZonoObsSet obs_set;
  for (int col = 0; col < num_cols; ++col) {
    const double x0 = data[0][col];
    const double y0 = data[1][col];
    const double h0 = data[2][col];
    const double vel = data[3][col];
    const double len = data[4][col];
    const double width = data[5][col];
    obs_set.PushObs(::roahm::DynObs{x0, y0, h0, vel, len, width, obs_type});
  }
  return obs_set;
}

std::vector<::roahm::MaybeRisky>
ProcessChooseable(const matlab::data::Array& arr) {
  fmt::print("Chooseable Arr Num Elts: {}\n", arr.getNumberOfElements());
  if (arr.isEmpty() or arr.getNumberOfElements() == 0) {
    return {};
  }

  if (arr.getType() != matlab::data::ArrayType::STRUCT) {
    fmt::print("CHOOSEABLE ERR 0\n");
    throw std::runtime_error("Expected struct array for maybe risky obs");
  }

  constexpr std::string_view kMaybeRiskyDynObsFieldName{"dyn_obs"};
  constexpr std::string_view kMaybeRiskyMuSigmaFieldName{"mu_sigma"};
  constexpr int kNumExpectedFields = 2;
  constexpr std::array<std::string_view, kNumExpectedFields>
      kExpectedFieldNames{kMaybeRiskyDynObsFieldName,
                          kMaybeRiskyMuSigmaFieldName};
  const matlab::data::StructArray& struct_arr{arr};
  if (struct_arr.getNumberOfFields() != kNumExpectedFields) {
    fmt::print("CHOOSEABLE ERR 1\n");
    throw std::runtime_error(fmt::format("Expected {} fields but got {}",
                                         kNumExpectedFields,
                                         struct_arr.getNumberOfFields()));
  }

  {
    int num_found_fields{0};
    for (const auto& field_name_mat_str : struct_arr.getFieldNames()) {
      for (const auto& expected_field_name : kExpectedFieldNames) {
        if (std::string{field_name_mat_str} ==
            std::string{expected_field_name}) {
          ++num_found_fields;
          break;
        }
      }
    }

    if (num_found_fields != kNumExpectedFields) {
      fmt::print("CHOOSEABLE ERR 2\n");
      throw std::runtime_error(fmt::format("Found {} of {} expected fields",
                                           num_found_fields,
                                           kNumExpectedFields));
    }
  }

  std::vector<roahm::MaybeRisky> ret;
  ret.reserve(struct_arr.getNumberOfElements());
  fmt::print("Struct Arr Num Elts: {}\n", struct_arr.getNumberOfElements());
  for (const auto& maybe_risky_mat_struct : struct_arr) {
    const auto& dyn_obs_portion = ProcessDynObs(
        maybe_risky_mat_struct[std::string{kMaybeRiskyDynObsFieldName}],
        ::roahm::DynObs::ObstacleType::kDynamicObs);
    const auto& risk_portion = ProcessAlwaysRiskyInput(
        maybe_risky_mat_struct[std::string{kMaybeRiskyMuSigmaFieldName}]);
    if (dyn_obs_portion.GetNumObs() != 1) {
      fmt::print("CHOOSEABLE ERR 3\n");
      throw std::runtime_error(
          fmt::format("Got {} dyn obs in single maybe risky instead of 1",
                      dyn_obs_portion.GetNumObs()));
    }
    if (risk_portion.size() != 1) {
      fmt::print("CHOOSEABLE ERR 4\n");
      throw std::runtime_error(
          fmt::format("Got {} mu sigma sets in single maybe risky instead of 1",
                      risk_portion.size()));
    }
    ret.emplace_back(dyn_obs_portion.GetSingleDynObs(0), risk_portion.at(0));
  }
  fmt::print("Ret Num Chooseable: {}\n", ret.size());
  return ret;
}

inline bool IsArgumentString(matlab::mex::ArgumentList& args, const int idx) {
  return (args.size() > idx) and
         (args[idx].getType() == matlab::data::ArrayType::MATLAB_STRING) and
         (args[idx].getDimensions().size() == 2) and
         (args[idx].getDimensions()[0] == 1) and
         (args[idx].getDimensions()[1] == 1);
}
inline bool IsArgumentListOfSizeAtLeast(matlab::mex::ArgumentList& args,
                                        const int num_min_expected) {
  return args.size() >= num_min_expected;
}
inline void ValidateArgumentSizeOfAtLeast(matlab::mex::ArgumentList& args,
                                          const int num_min_expected) {
  if (not IsArgumentListOfSizeAtLeast(args, num_min_expected)) {
    throw std::runtime_error("Expected at least " +
                             std::to_string(num_min_expected) + " arguments!");  }
}
inline void ValidateArgumentIsSingleString(matlab::mex::ArgumentList& args,
                                           const int idx) {
  ValidateArgumentSizeOfAtLeast(args, idx);
  if (not IsArgumentString(args, idx)) {
    throw std::runtime_error("Expected single string in argument " +
                             std::to_string(idx) + "!");
  }
}

inline std::string ValidateAndLoadSingleString(matlab::mex::ArgumentList& args,
                                               const int idx) {
  ValidateArgumentSizeOfAtLeast(args, idx);
  ValidateArgumentIsSingleString(args, idx);
  const matlab::data::StringArray ai = args[idx];
  return std::string(args[idx][0]);
}
inline bool IsArgumentDoubleMat(matlab::mex::ArgumentList& args,
                                const int idx) {
  return (args.size() > idx) and
         (args[idx].getType() == matlab::data::ArrayType::DOUBLE);
}

inline void ValidateArgumentIsDoubleMat(matlab::mex::ArgumentList& args,
                                        const int idx) {
  ValidateArgumentSizeOfAtLeast(args, idx);
  if (not IsArgumentDoubleMat(args, idx)) {
    throw std::runtime_error(
        "Expected matrix of doubles for argument at index " +
        std::to_string(idx));
  }
}
inline bool DoubleMatHasSize(matlab::mex::ArgumentList& args, const int idx,
                             const int rows, const int cols) {
  return IsArgumentDoubleMat(args, idx) and
         (args[idx].getDimensions().size() == 2) and
         (args[idx].getDimensions()[0] == rows) and
         (args[idx].getDimensions()[1] == cols);
}


void ValidateArgumentIsDoubleMatWithSize(matlab::mex::ArgumentList& args,
                                         const int idx, const int rows,
                                         const int cols) {
  ValidateArgumentSizeOfAtLeast(args, idx);
  ValidateArgumentIsDoubleMat(args, idx);
  if (not DoubleMatHasSize(args, idx, rows, cols)) {
    if (args[idx].getDimensions().size() != 2) {
      throw std::runtime_error(
          "Expected 2 dimensions, but got " +
          std::to_string(args[idx].getDimensions().size()));
    }
    throw std::runtime_error(
        "Expected [" + std::to_string(rows) + " x " + std::to_string(cols) +
        "] double matrix but got [" +
        std::to_string(args[idx].getDimensions()[0]) + " x " +
        std::to_string(args[idx].getDimensions()[1]) + "] instead");
  }
}
inline double ValidateAndLoadArgumentSingleDoubleValue(
    matlab::mex::ArgumentList& args, const int idx) {
  ValidateArgumentIsDoubleMatWithSize(args, idx, 1, 1);
  return double{args[idx][0][0]};
}

::roahm::RiskProblemDescription
ValidateAndLoadInputs(matlab::mex::ArgumentList& outputs,
                      matlab::mex::ArgumentList& inputs) {
  // Expectd Inputs TODO should this be a struct:
  // 1. (Predicted) State (x, y, theta, u, v, r)
  // 2. A waypoint TODO accept multiple?
  // 3. chooseable obs {mu_sigma_10ms, mu_sigma_20ms, dyn_obs}

  // 4. TODO real obstacles

  // Expected Outputs:
  // 1. Array TODO
  constexpr int kNumExpectedInputs = 7;
  constexpr int kNumExpectedOutputs = 1;

  constexpr int kStateVarInputNum = 0;
  constexpr int kWaypointInputNum = 1;
  constexpr int kChooseableObsInputNum = 2;
  constexpr int kAlwaysNonRiskyObsInputNum = 3;
  constexpr int kAlwaysRiskyObsInputNum = 4;

  static_assert(kNumExpectedInputs > kStateVarInputNum,
                "Not enough inputs expected to access rover state variable");
  static_assert(kNumExpectedInputs > kWaypointInputNum,
                "Not enough inputs expected to access waypoint variable");
  static_assert(
      kNumExpectedInputs > kChooseableObsInputNum,
      "Not enough inputs expected to access chooseable obs variable variable");
  static_assert(kNumExpectedInputs > kAlwaysNonRiskyObsInputNum,
                "Not enough inputs expected to always non risky obs variable");
  static_assert(kNumExpectedInputs > kAlwaysRiskyObsInputNum,
                "Not enough inputs expected to always risky obs variable");

  if (inputs.size() != kNumExpectedInputs) {
    throw std::runtime_error(fmt::format("Expected {} inputs, but got {}",
                                         kNumExpectedInputs, inputs.size()));
  }
  if (outputs.size() != kNumExpectedOutputs) {
    throw std::runtime_error(fmt::format("Expected {} outputs, but got {}",
                                         kNumExpectedOutputs, outputs.size()));
  }

  ::roahm::RiskProblemDescription ret{};
  fmt::print("Loading Rover State...\n");
  ret.rover_state_ = LoadRoverStateFromMatlab(inputs[kStateVarInputNum]);
  fmt::print("Loading Waypoint...\n");
  ret.waypoint_global_no_mirror_ =
      LoadWaypointFromMatlab(inputs[kWaypointInputNum]);
  fmt::print("Loading Always Risky...\n");
  ret.always_risky_obs_ =
      ProcessAlwaysRiskyInput(inputs[kAlwaysRiskyObsInputNum]);
  fmt::print("Loading Always Non Risky...\n");
  ret.always_non_risky_obs_ =
      ProcessDynObs(inputs[kAlwaysNonRiskyObsInputNum],
                    roahm::DynObs::ObstacleType::kStaticBoundary);
  fmt::print("Loaded {} Always Non Risky\n",
             ret.always_non_risky_obs_.GetNumObs());
  fmt::print("Loading Maybe Risky...\n");
  ret.maybe_risky_obs_ = ProcessChooseable(inputs[kChooseableObsInputNum]);
  fmt::print("Returning...\n");
  return ret;
}

class MexFunction : public matlab::mex::Function {
private:
  std::string last_frs_fpath_;
  ::roahm::RiskRtd risk_rtd_;

public:
  MexFunction() : last_frs_fpath_{"/data/cpp_processed_CUDA_FRS_30-Aug-2022_lessFRS_3rnd_grid24_t_fix.DBG"}, risk_rtd_{last_frs_fpath_} {}
  void operator()(matlab::mex::ArgumentList outputs,
                  matlab::mex::ArgumentList inputs) {
    using ::roahm::GetDeltaS;
    using ::roahm::Tick;
    try {
      constexpr bool kPlotThings{false};
      matlab::data::ArrayFactory factory;
      matlab::data::TypedArray<double> arr_out =
          factory.createArray<double>({8});

      fmt::print("[RRM] Getting Engine\n");
      std::shared_ptr<matlab::engine::MATLABEngine> engine_ptr = getEngine();
      // Setup if we want to do any plotting
      if (kPlotThings) {
        roahm::plot_utils::SetupEngine(engine_ptr);
      }

      fmt::print("[RRM] Getting Inputs\n");
      auto risk_inputs{ValidateAndLoadInputs(outputs, inputs)};
      const std::string new_frs_fpath{ValidateAndLoadSingleString(inputs, 5)};
      const double risk_threshold{ValidateAndLoadArgumentSingleDoubleValue(inputs, 6)};
      fmt::print("[LOADED RISK THRESHOLD] {}\n", risk_threshold);
      if (new_frs_fpath != last_frs_fpath_) {
        fmt::print("[SETTING NEW FRS FPATH] {}\n", new_frs_fpath);
	last_frs_fpath_ = new_frs_fpath;
        fmt::print("[LOADING FRS @ FPATH] {}\n", new_frs_fpath);
	risk_rtd_ = ::roahm::RiskRtd{last_frs_fpath_};
        fmt::print("[LOADED FRS @ FPATH] {}\n", new_frs_fpath);
      }
      fmt::print("State: {} {}\n", risk_inputs.rover_state_.x_,
                 risk_inputs.rover_state_.y_);
      fmt::print("[RRM] Getting Outputs\n");
      // fmt::print("[MEX] Ran Planning Iteration\n");
      // fmt::print("[MEX] RISK RTD Took: {}\n",
      // risk_outputs.total_time_seconds_);

      fmt::print("[RRM] Configuring Slicing\n");
      const ::roahm::plot_utils::PlotInfo plot_info_sliced{
          ::roahm::plot_utils::PlotColor{0.01, 0.0, 1.0, 0.0}, 0.5, false,
          ::roahm::plot_utils::PlotLineType::Continuous()};
      fmt::print("[RRM] Loading FRS\n");
      const auto frs_full{roahm::LoadFrsBinary(last_frs_fpath_)};
      //const auto frs_full{roahm::LoadFrsBinary(
      //    "/data/"
      //    "cpp_processed_CUDA_FRS_30-Aug-2022_lessFRS_3rnd_grid24_t_fix.DBG")};
      fmt::print("[RRM] Loaded FRS\n");
      // roahm::plot_impls::PlotFlZonoObsSet(risk_inputs.always_non_risky_obs_);
      fmt::print("Num Non Risky:   {}\n",
                 risk_inputs.always_non_risky_obs_.GetNumObs());
      fmt::print("Num Risky:       {}\n", risk_inputs.always_risky_obs_.size());
      fmt::print("Num Maybe Risky: {}\n", risk_inputs.maybe_risky_obs_.size());
      roahm::plot_utils::PlotInfo plot_info{
          roahm::plot_utils::PlotColor{0.0, 0.0, 0.0, 0.0}, 1, true,
          roahm::plot_utils::PlotLineType::Continuous(), std::nullopt};
      if (kPlotThings) {
        // roahm::PlotVehrs(frs_full.megas_.at(0).dir_.at(0).at(0),
        //                  risk_inputs.rover_state_.GetXYH(), false,
        //                  plot_info_sliced);

        // roahm::plot_utils::TurnHoldOn();
        roahm::plot_utils::OpenFigure(34);
        roahm::plot_utils::ClearAxes();
        // roahm::plot_utils::TurnHoldOn();
        // roahm::plot_utils::OpenFigure(13);
        roahm::plot_utils::TurnHoldOn();
        roahm::plot_impls::PlotPointXYH(risk_inputs.rover_state_.GetXYH(), 4.0,
                                        1.0, plot_info);
        roahm::plot_impls::PlotPointXYH(risk_inputs.waypoint_global_no_mirror_,
                                        4.0, 1.0, plot_info);
        for (const auto& maybe_obs : risk_inputs.maybe_risky_obs_) {
          roahm::plot_utils::TurnHoldOn();
          roahm::plot_impls::PlotFlZonoObs(maybe_obs.dyn_obs_);
        }
        roahm::plot_utils::SetAxisEqual();
      }

      fmt::print("[RMM] Pre Planning ITER\n");
      const auto t0 = Tick();
      std::vector<roahm::RiskRtdPlanningOutputsIpopt> risk_outputs_original;
      // CHALLEN_NOTE abl_meth is 1-indexed to the list of methods 
      // you posted on Slack.
      const int abl_meth = -1;
      if (abl_meth == 1) {
        const std::uint32_t randomization_seed{0};
        const ::roahm::TimeMicroseconds dt{
            ::roahm::TimeMicroseconds::Get10ms()};
        //const double risk_threshold{0.05};
        const std::int64_t num_trajectory_samples{1000};
        const roahm::risk_comparisons::EnvFootprints footprints{
            roahm::risk_comparisons::EnvFootprints::Default()};
        const int num_param_samples{1};
        fmt::print("ABL_METH: 1\n");
        const auto sampled_outputs =
            roahm::risk_comparisons::MonteCarloDiscreteSampling27(
                frs_full, risk_inputs, randomization_seed,
                num_trajectory_samples, dt, footprints, risk_threshold, num_param_samples);
        fmt::print("OUT: {}\n",
                   sampled_outputs.at(0).final_cost_.value_or(-1.00));
        std::vector<roahm::RiskRtdPlanningOutputsIpopt> fake_ipopt;
        for (const auto& samp_out : sampled_outputs) {
          roahm::RiskRtdPlanningOutputsIpopt nsamp{
              samp_out, Ipopt::ApplicationReturnStatus::Internal_Error,
              std::nullopt};
          fake_ipopt.push_back(nsamp);
        }
        risk_outputs_original = fake_ipopt;
      } else if (abl_meth == 3) {
        const std::uint32_t randomization_seed{0};
        const ::roahm::TimeMicroseconds dt{
            ::roahm::TimeMicroseconds::Get10ms()};
        //const double risk_threshold{0.05};
        const std::int64_t num_trajectory_samples{1000};
        const roahm::risk_comparisons::EnvFootprints footprints{
            roahm::risk_comparisons::EnvFootprints::Default()};
        const auto sampled_outputs =
            roahm::risk_comparisons::MonteCarloDiscreteSamplingOptE(
                frs_full, risk_inputs, randomization_seed,
                num_trajectory_samples, dt, footprints, risk_threshold);
        fmt::print("OUT: {}\n",
                   sampled_outputs.at(0).final_cost_.value_or(-1.00));
        std::vector<roahm::RiskRtdPlanningOutputsIpopt> fake_ipopt;
        for (const auto& samp_out : sampled_outputs) {
          roahm::RiskRtdPlanningOutputsIpopt nsamp{
              samp_out, Ipopt::ApplicationReturnStatus::Internal_Error,
              std::nullopt};
          fake_ipopt.push_back(nsamp);
        }
        risk_outputs_original = fake_ipopt;
      } else if (abl_meth == 2) {
        risk_outputs_original = risk_rtd_.RunPlanningIteration<::roahm::monte_carlo_27_ipopt::MonteCarlo27IpoptProblem>(risk_inputs, true, risk_threshold);
      } else if (abl_meth == 4) {
      risk_outputs_original = risk_rtd_.RunPlanningIteration<::roahm::monte_carlo_ipopt::MonteCarloIpoptProblem>(risk_inputs, true, risk_threshold);
      } else {
        risk_outputs_original = risk_rtd_.RunPlanningIteration<::roahm::risk_rtd_ipopt_problem::RiskRtdIpoptProblem>(risk_inputs, true, risk_threshold);
      }
      auto risk_outputs = risk_outputs_original;
      const auto t1 = Tick();

      // std::vector<roahm::plot_utils::PlotColor> color_its{
      //     roahm::plot_utils::PlotColor{0.5, 0.0, 0.0, 0.0},
      //     roahm::plot_utils::PlotColor{0.5, 1.0, 0.0, 0.0},
      //     roahm::plot_utils::PlotColor{0.5, 0.0, 1.0, 0.0},
      //     roahm::plot_utils::PlotColor{0.5, 0.0, 0.0, 1.0},
      //     roahm::plot_utils::PlotColor{0.5, 1.0, 1.0, 0.0},
      //     roahm::plot_utils::PlotColor{0.5, 0.0, 1.0, 1.0},
      //     roahm::plot_utils::PlotColor{0.5, 1.0, 0.0, 1.0},
      // };

      fmt::print("[RRM] Loading FRS\n");
      const double delta_y_wp{risk_inputs.waypoint_global_no_mirror_.y_ -
                              risk_inputs.rover_state_.y_};
      const bool is_waypoint_in_same_lane{std::abs(delta_y_wp) <= 0.1};
      {
        int color_idx{0};
        int idx{0};
        int success_count{0};
        fmt::print("[RRM] Iterating\n");
        for (const auto& full_opt_inf : risk_outputs) {
          ++idx;
          success_count += full_opt_inf.found_feasible_;
          // const ::roahm::plot_utils::PlotInfo failed_plot_info{
          //     ::roahm::plot_utils::PlotColor{0.5, 1.0, 0.0, 0.0}, 0.05,
          //     false,
          //     ::roahm::plot_utils::PlotLineType::Continuous(), std::nullopt};
          // const ::roahm::plot_utils::PlotInfo success_plot_info{
          //     ::roahm::plot_utils::PlotColor{0.5, 0.0, 1.0, 0.0}, 0.05,
          //     false,
          //     ::roahm::plot_utils::PlotLineType::Continuous(), std::nullopt};
          // auto plot_info = full_opt_inf.found_feasible_ ? success_plot_info
          //                                               : failed_plot_info;
          const auto& curr_frs = frs_full.GetVehrs(full_opt_inf.sel_inf_);
          if (true or (not full_opt_inf.found_feasible_)) {
            fmt::print(R"""([ RES ] Problem {} / {} (1-idx)
  Mirror:        {}
  ManuType:      {}
  Feasible:      {}
  Solve Success: {}  
  Status:        {}
  Param:         {} in [{}, {}]
  Cost:          {}
  Solver Status: {}
)""",
                       idx, risk_outputs.size(), full_opt_inf.sel_inf_.mirror_,
                       ToString(full_opt_inf.sel_inf_.manu_type_),
                       full_opt_inf.found_feasible_,
                       full_opt_inf.solve_succeeded_,
                       ::roahm::ToString(full_opt_inf.ipopt_status_),
                       ::roahm::StdToStringOpt(full_opt_inf.final_param_),
                       curr_frs.GetTrajParamMin(), curr_frs.GetTrajParamMax(),
                       ::roahm::StdToStringOpt(full_opt_inf.final_cost_),
                       ::roahm::ToStringOpt<Ipopt::SolverReturn>(
                           full_opt_inf.final_solver_status_));

            fmt::print("Center val: {}\n", curr_frs.GetCenterK());
            fmt::print("Idx: {} Color {}\n", idx, color_idx);
            // const auto new_color =
            //     color_its.at((++color_idx) % color_its.size());
            //  plot_info.color_ = new_color;
            if (full_opt_inf.found_feasible_) {
              plot_info.color_ = {0.5, 0.0, 1.0, 0.0};
            } else {
              plot_info.color_ = {0.3, 1.0, 0.0, 0.0};
            }
            const double param_val_unmirrored =
                full_opt_inf.final_param_.value_or(curr_frs.GetTrajParamMin());
            if (kPlotThings and IsSpd(curr_frs.GetManuType()) and
                not full_opt_inf.found_feasible_ and
                not full_opt_inf.sel_inf_.mirror_) {
              fmt::print("Is Spd\n");
              // PlotVehrs(
              //     curr_frs.SliceAtParam(risk_inputs.rover_state_.u_,
              //                           risk_inputs.rover_state_.v_,
              //                           risk_inputs.rover_state_.r_,
              //                           full_opt_inf.final_param_.value_or(
              //                               curr_frs.GetTrajParamMin())),
              //     risk_inputs.rover_state_.GetXYH(),
              //     full_opt_inf.sel_inf_.mirror_, plot_info);
              ::roahm::plot_utils::TurnHoldOn();
              // PlotVehrs(
              //     curr_frs.SliceAtParam(risk_inputs.rover_state_.u_,
              //                           risk_inputs.rover_state_.v_,
              //                           risk_inputs.rover_state_.r_,
              //                           full_opt_inf.final_param_.value_or(
              //                               curr_frs.GetTrajParamMax())),
              //     risk_inputs.rover_state_.GetXYH(),
              //     full_opt_inf.sel_inf_.mirror_, plot_info);
              PlotVehrs(
                  curr_frs.SliceAtParam(risk_inputs.rover_state_.u_,
                                        risk_inputs.rover_state_.v_,
                                        risk_inputs.rover_state_.r_,
                                        full_opt_inf.final_param_.value_or(
                                            curr_frs.GetCenterK())),
                  risk_inputs.rover_state_.GetXYH(),
                  full_opt_inf.sel_inf_.mirror_, plot_info);
            }
          }
        }
        fmt::print("{} successful / {}\n", success_count, risk_outputs.size());
        fmt::print("waypoint_in_same_lane: {}\n", is_waypoint_in_same_lane);
        fmt::print("global waypoint: {} {} {}\n",
                   risk_inputs.waypoint_global_no_mirror_.x_,
                   risk_inputs.waypoint_global_no_mirror_.y_,
                   risk_inputs.waypoint_global_no_mirror_.h_);
        const auto local_tmp =
            risk_inputs.waypoint_global_no_mirror_.RelativeToEgoFrameNoMirror(
                risk_inputs.rover_state_.GetXYH());
        fmt::print("local tmp: {} {} {}\n", local_tmp.x_, local_tmp.y_,
                   local_tmp.h_);
        fmt::print("start point: {} {} {}\n", risk_inputs.rover_state_.x_,
                   risk_inputs.rover_state_.y_,
                   risk_inputs.rover_state_.GetHeading());

        if (success_count > 0) {
          fmt::print("Success count > 0\n");
          // Found feasible
          const double u0{risk_inputs.rover_state_.u_};
          const double h0{risk_inputs.rover_state_.GetHeading()};
          const auto rtd_out_comparison_operator =
              [u0, h0, delta_y_wp, is_waypoint_in_same_lane](
                  const ::roahm::RiskRtdPlanningOutputs& p0,
                  const ::roahm::RiskRtdPlanningOutputs& p1) -> bool {
            // p0 < p1  ===>  p0 has higher priority
            // (true)
            constexpr bool kPrioritizeP0{true};
            constexpr bool kPrioritizeP1{false};
            static_assert(kPrioritizeP0 != kPrioritizeP1);

            /*
            if (p0.found_feasible_ and (not p1.found_feasible_)) {
              return kPrioritizeP0;
            } else if (p1.found_feasible_ and (not p0.found_feasible_)) {
              return kPrioritizeP1;
            } else if ((not p0.found_feasible_) and (not p1.found_feasible_)) {
              // TODO is this a total ordering?
              const bool idx_tiebreaker =
                  (p0.sel_inf_.idxu0_ < p1.sel_inf_.idxu0_) or
                  (p0.sel_inf_.idx0_ < p1.sel_inf_.idx0_) or
                  (p0.sel_inf_.idx1_ < p1.sel_inf_.idx1_);
              return idx_tiebreaker;
            }
            */

            {
              const bool p0_has_param{p0.final_param_.has_value()};
              const bool p1_has_param{p1.final_param_.has_value()};
              const bool p0_has_cost{p0.final_cost_.has_value()};
              const bool p1_has_cost{p1.final_cost_.has_value()};
              const bool p0_has_location{p0.final_location_.has_value()};
              const bool p1_has_location{p1.final_location_.has_value()};
              if ((not p0_has_param) or (not p1_has_param)) {
                throw std::runtime_error("Feasible with no params");
              }
              if ((not p0_has_cost) or (not p1_has_cost)) {
                throw std::runtime_error("Feasible with no cost");
              }
              if ((not p0_has_location) or (not p1_has_location)) {
                throw std::runtime_error("Feasible with no location");
              }
            }
            const double p0_final_param{p0.final_param_.value()};
            const double p1_final_param{p1.final_param_.value()};
            const double p0_cost{p0.final_cost_.value()};
            const double p1_cost{p1.final_cost_.value()};
            const ::roahm::PointXYH p0_location_local{
                p0.final_location_.value()};
            const ::roahm::PointXYH p1_location_local{
                p1.final_location_.value()};
            const double p0_mirror_mult{p0.sel_inf_.mirror_ ? -1.0 : 1.0};
            const double p1_mirror_mult{p1.sel_inf_.mirror_ ? -1.0 : 1.0};
            const double p0_delta_y{p0_location_local.y_ * p0_mirror_mult};
            const double p1_delta_y{p1_location_local.y_ * p1_mirror_mult};
            // const double p0_abs_delta_y{std::abs(p0_dy)};
            // const double p1_abs_delta_y{std::abs(p1_dy)};
            const double p0_signed_dist_to_wp_y{delta_y_wp - p0_delta_y};
            const double p1_signed_dist_to_wp_y{delta_y_wp - p1_delta_y};
            const double p0_abs_dist_to_wp_y{std::abs(p0_signed_dist_to_wp_y)};
            const double p1_abs_dist_to_wp_y{std::abs(p1_signed_dist_to_wp_y)};

            // TODO check heading

            constexpr bool kIsRover{true};
            constexpr double kWantToSpeedUpThresholdRover{0.0};
            constexpr double kWantToSpeedUpThresholdCar{10.6};
            constexpr double kWantToSpeedUpThreshold{kIsRover ? kWantToSpeedUpThresholdRover : kWantToSpeedUpThresholdCar}; // [m/s]
            const bool really_want_to_speed_up{u0 < kWantToSpeedUpThreshold};

            const bool is_p0_lan{IsLan(p0.sel_inf_.manu_type_)};
            const bool is_p1_lan{IsLan(p1.sel_inf_.manu_type_)};
            const bool is_p0_spd{IsSpd(p0.sel_inf_.manu_type_)};
            const bool is_p1_spd{IsSpd(p1.sel_inf_.manu_type_)};
            const bool is_p0_dir{IsDir(p0.sel_inf_.manu_type_)};
            const bool is_p1_dir{IsDir(p1.sel_inf_.manu_type_)};

            const bool is_p0_spd_up{(is_p0_spd) and (p0_final_param > u0)};
            const bool is_p1_spd_up{(is_p1_spd) and (p1_final_param > u0)};

            const bool is_p0_slow_down{(is_p0_spd) and (p0_final_param < u0)};
            const bool is_p1_slow_down{(is_p1_spd) and (p1_final_param < u0)};

            if (std::abs(h0) > 0.1) {
              const double delta_h_p0{(p0.sel_inf_.mirror_ ? -1.0 : 1.0) *
                                      p0.final_location_.value().h_};
              const double delta_h_p1{(p1.sel_inf_.mirror_ ? -1.0 : 1.0) *
                                      p1.final_location_.value().h_};
              const double h_p0_global{h0 + delta_h_p0};
              const double h_p1_global{h0 + delta_h_p1};
              const double abs_h_p0_global{std::abs(h_p0_global)};
              const double abs_h_p1_global{std::abs(h_p1_global)};
              const bool delta_h_p0_improves{abs_h_p0_global < std::abs(h0)};
              const bool delta_h_p1_improves{abs_h_p1_global < std::abs(h0)};
              if (delta_h_p0_improves and delta_h_p1_improves) {
                if (abs_h_p0_global < abs_h_p1_global) {
                  return kPrioritizeP0;
                } else if (abs_h_p1_global < abs_h_p0_global) {
                  return kPrioritizeP1;
                }
              } else if (delta_h_p0_improves and (not delta_h_p1_improves)) {
                return kPrioritizeP0;
              } else if (delta_h_p1_improves and (not delta_h_p0_improves)) {
                return kPrioritizeP1;
              }
            }

            // If we are going slow, prioritize a speed change
            if (really_want_to_speed_up or is_waypoint_in_same_lane) {
              if (is_p0_spd_up and is_p1_spd_up) {
                // If both are speed up, check which one speeds up more
                if (p0_cost < p1_cost) {
                  return kPrioritizeP0;
                } else if (p1_cost < p0_cost) {
                  return kPrioritizeP1;
                }
              } else if (is_p0_spd_up) {
                // If p0 speeds up, but p1 does not, prioritize p0
                return kPrioritizeP0;
              } else if (is_p1_spd_up) {
                // If p1 speeds up, but p0 does not, prioritize p1
                return kPrioritizeP1;
              } else if (is_p0_slow_down) {
                // If p0 slows down, but p1 maintains, prioritize p1
                return kPrioritizeP1;
              } else if (is_p1_slow_down) {
                // If p1 slows down, but p0 maintains, prioritize p0
                return kPrioritizeP0;
              }
              // If we get here, neither are speed up or slow down
            }

            if (is_waypoint_in_same_lane) {
              if (p0_abs_dist_to_wp_y < p1_abs_dist_to_wp_y) {
                return kPrioritizeP0;
              } else if (p1_abs_dist_to_wp_y < p0_abs_dist_to_wp_y) {
                return kPrioritizeP1;
              }
            }

            if (not is_waypoint_in_same_lane) {
              const bool is_p0_ok{is_p0_lan or is_p0_spd};
              const bool is_p1_ok{is_p1_lan or is_p1_spd};
              if (is_p0_ok and is_p1_ok) {
                if (p0_abs_dist_to_wp_y < p1_abs_dist_to_wp_y) {
                  return kPrioritizeP0;
                } else if (p1_abs_dist_to_wp_y < p0_abs_dist_to_wp_y) {
                  return kPrioritizeP1;
                }
              } else if (is_p0_ok) {
                return kPrioritizeP0;
              } else if (is_p1_ok) {
                return kPrioritizeP1;
              }
            }

            // Left the full if statement here for clarity, and in case we
            // change sort order (and kPrioritizeP0, kPrioritizeP1)
            if (p0_cost < p1_cost) {
              return kPrioritizeP0;
            }
            return kPrioritizeP1;

            // Prioritize Lane Changes
            // const bool is_p0_more_sketchy = (p0_in_sketchy_lane and (not
            // p1_in_sketchy_lane)); const bool is_p1_more_sketchy =
            // (p1_in_sketchy_lane and (not p0_in_sketchy_lane));
            /*
            const bool is_p0_significant_lan{
                is_p0_lan and
                (std::abs(p0_final_param) >= kSignificantLaneChangeMinAbs)};
            const bool is_p1_significant_lan{
                is_p1_lan and
                (std::abs(p1_final_param) >= kSignificantLaneChangeMinAbs)};

            if (not is_waypoint_in_same_lane) {
              if (is_p0_significant_lan and (not is_p1_significant_lan)) {
                return kPrioritizeP0;
              } else if (is_p1_significant_lan and
                         (not is_p0_significant_lan)) {
                return kPrioritizeP1;
              } else if (is_p0_lan and (not is_p1_lan)) {
                return kPrioritizeP0;
              } else if (is_p1_lan and (not is_p0_lan)) {
                return kPrioritizeP1;
              }
            }
            */

            // Try not to actively slow down, in general
            /*
            if (is_p1_slow_down and (not is_p0_slow_down)) {
              return kPrioritizeP0;
            } else if (is_p0_slow_down and (not is_p1_slow_down)) {
              return kPrioritizeP1;
            }
            */

            /*
            if (is_p0_lan and (not is_p1_lan)) {
              return kPrioritizeP0;
            } else if (is_p1_lan and (not is_p0_lan)) {
              return kPrioritizeP1;
            }
            */

            /*
            if (is_p1_dir and (not is_p0_dir)) {
              return kPrioritizeP0;
            } else {
              return kPrioritizeP1;
            }

            return p0_cost < p1_cost;
            */
          };

          for (int i = 0; i < risk_outputs.size(); ++i) {
            if (not risk_outputs.at(i).found_feasible_) {
              risk_outputs.erase(risk_outputs.begin() + i);
              --i;
            }
          }

          std::sort(risk_outputs.begin(), risk_outputs.end(),
                    rtd_out_comparison_operator);
          if (risk_outputs.size() > 1) {
            std::sort(risk_outputs.begin(), risk_outputs.begin() + 2,
                      rtd_out_comparison_operator);
          }

          for (int i = 0; i < risk_outputs.size(); ++i) {
            const auto& p0 = risk_outputs.at(i);
            fmt::print(
                "[{}]\n Type: {}\n P: {} of [{}, {}]\n WPLMY: {}\n FLY: {}\n",
                i, ToString(p0.sel_inf_.manu_type_), p0.final_param_.value(),
                frs_full.GetVehrs(p0.sel_inf_).GetTrajParamMin(),
                frs_full.GetVehrs(p0.sel_inf_).GetTrajParamMax(),
                p0.waypoint_local_frame_.y_, p0.final_location_.value().y_);
          }
          if (not risk_outputs.front().found_feasible_) {
            throw std::runtime_error(
                "Sorted, success count > 0, but front doesn't have feasible?");
          }
          fmt::print("Constructed\n");
          /*
          for (int i = 0; i < risk_outputs.size(); ++i) {
            const auto& curr_out = risk_outputs.at(i);
            if (curr_out.found_feasible_) {
              if (curr_out.final_cost_.has_value()) {
                constexpr double kInf = std::numeric_limits<double>::infinity();
                const double curr_cost = curr_out.final_cost_.value();
                const bool have_prev_min = not(min_cost.has_value());
                const bool curr_lower_cost =
                    (curr_cost < min_cost.value_or(kInf));
                const bool curr_is_lan =
                    curr_out.sel_inf_.manu_type_ == roahm::ManuType::kLanChange;
                const bool prev_is_lan =
                    min_manu_type.value_or(::roahm::ManuType::kNone) ==
                    ::roahm::ManuType::kLanChange;
                if (have_prev_min or curr_lower_cost or
                    (curr_is_lan and (not prev_is_lan))) {
                  min_cost = curr_cost;
                  min_param_no_mirror = curr_out.final_param_;
                  min_manu_type = curr_out.sel_inf_.manu_type_;
                  min_mirror = curr_out.sel_inf_.mirror_;
                  min_u0 = curr_out.u0_;
                  min_sel_info = curr_out.sel_inf_;
                }
              }
            }
          }
          */

          const auto& min_elt = risk_outputs.front();
          std::optional<double> min_cost = min_elt.final_cost_;
          std::optional<double> min_param_no_mirror = min_elt.final_param_;
          std::optional<double> min_u0 = min_elt.u0_;
          std::optional<::roahm::ManuType> min_manu_type =
              min_elt.sel_inf_.manu_type_;
          std::optional<::roahm::FrsSelectInfo> min_sel_info = min_elt.sel_inf_;
          const bool min_mirror{min_elt.sel_inf_.mirror_};
          int min_manu_type_int = 0;
          if (min_manu_type == ::roahm::ManuType::kSpdChange) {
            min_manu_type_int = 1;
          } else if (min_manu_type == ::roahm::ManuType::kDirChange) {
            min_manu_type_int = 2;
          } else if (min_manu_type == ::roahm::ManuType::kLanChange) {
            min_manu_type_int = 3;
          }
          fmt::print("Checkpoint 00\n");
          constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
          arr_out[0] = true;
          arr_out[1] = min_param_no_mirror.value_or(kNaN);
          arr_out[2] = min_mirror;
          arr_out[3] = min_manu_type_int;
          arr_out[4] = min_u0.value_or(kNaN);
          fmt::print("Checkpoint 01\n");
          if (not min_sel_info.has_value()) {
            throw std::runtime_error("No sel info on final value");
          }
          // const auto sel_info = min_sel_info.value_or(
          // roahm::FrsSelectInfo{roahm::ManuType::kNone, 0, 0, 0, false});
          const auto sel_info = min_sel_info.value();
          arr_out[5] = sel_info.idxu0_;
          arr_out[6] = sel_info.idx0_;
          arr_out[7] = sel_info.idx1_;
          if (min_param_no_mirror.has_value()) {
            const auto state_xyh{risk_inputs.rover_state_.GetXYH()};
            const auto min_param_no_mirror_val{min_param_no_mirror.value()};
            const auto got_vehrs{frs_full.GetVehrs(sel_info)};
            const auto sliced_vehrs{got_vehrs.SliceAtParam(
                risk_inputs.rover_state_.u_, risk_inputs.rover_state_.v_,
                risk_inputs.rover_state_.r_, min_param_no_mirror_val)};
            if (kPlotThings) {
              roahm::PlotVehrs(sliced_vehrs, state_xyh, min_mirror,
                               std::nullopt);
            }
          }
          if (::roahm::IsSpd(
                  min_manu_type.value_or(::roahm::ManuType::kNone))) {
            fmt::print("[RRM] Got Outputs [succ: {}] in {}sec, spd {}mps\n",
                       true, GetDeltaS(t1, t0),
                       min_param_no_mirror.value_or(kNaN));
          } else {
            fmt::print("[RRM] Got Outputs [succ: {}] in {}sec, spd {}mps\n",
                       true, GetDeltaS(t1, t0), risk_inputs.rover_state_.u_);
          }
        } else {
          fmt::print("No Successes, zeroing outputs\n");
          arr_out[0] = false;
          arr_out[1] = 0;
          arr_out[2] = 0;
          arr_out[3] = 0;
          arr_out[4] = 0;
          arr_out[5] = 0;
          arr_out[6] = 0;
          arr_out[7] = 0;
          fmt::print("[RRM] Got Outputs [succ: {}] in {}sec, spd {}mps\n",
                     false, GetDeltaS(t1, t0), risk_inputs.rover_state_.u_);
        }
      }
      fmt::print("Setting outputs[0]\n");
      outputs[0] = arr_out;
    } catch (std::runtime_error e) {
      fmt::print(stderr, "Runtime Error Caught: {}\n", e.what());
    } catch (std::out_of_range e) {
      fmt::print(stderr, "Out of Range Error Caught: {}\n", e.what());
    } catch (std::exception e) {
      fmt::print(stderr, "Exception Caught: {}\n", e.what());
    } catch (...) {
      fmt::print(stderr,
                 "Caught exception not inheriting from std::exception\n");
    }
  }
};
