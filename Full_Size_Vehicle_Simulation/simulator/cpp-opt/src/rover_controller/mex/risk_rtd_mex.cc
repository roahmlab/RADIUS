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
#include "dyn_obs.hpp"
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
#include "point_xyh.hpp"
#include "risk_rtd.hpp"
#include "risk_rtd_ipopt_problem.hpp"
#include "ros/time.h"
#include "rover_state.hpp"
#include "timing_util.hpp"
#include "waypoint.hpp"

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
  }

  constexpr int kNumExpectedDims = 2;
  const int num_dims = arr.getDimensions().size();
  //fmt::print("Number of dimension: {}\n", num_dims);
  if (num_dims != kNumExpectedDims) {
    throw std::runtime_error(fmt::format("Expected {} dimensions but got {}",
                                         kNumExpectedDims, num_dims));
  }

  const int num_mu_sigma_sets = arr.getDimensions()[0];
  //fmt::print("Dim[0]: {}\n", num_mu_sigma_sets);
  if (num_mu_sigma_sets <= 0) {
    return {};
  }

  constexpr int kNumExpectedMuSigmaDurations = 16;
  const int num_mu_sigma_durations = arr.getDimensions()[1];
  //fmt::print("Dim[1]: {}\n", num_mu_sigma_durations);
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
inline bool IsArgumentBoolMat(matlab::mex::ArgumentList& args,
                                const int idx) {
  return (args.size() > idx) and
         (args[idx].getType() == matlab::data::ArrayType::LOGICAL);
}
inline bool IsArgumentDoubleMat(matlab::mex::ArgumentList& args,
                                const int idx) {
  return (args.size() > idx) and
         (args[idx].getType() == matlab::data::ArrayType::DOUBLE);
}

inline void ValidateArgumentIsBoolMat(matlab::mex::ArgumentList& args,
                                        const int idx) {
  ValidateArgumentSizeOfAtLeast(args, idx);
  if (not IsArgumentBoolMat(args, idx)) {
    throw std::runtime_error(
        "Expected matrix of bools for argument at index " +
        std::to_string(idx));
  }
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
inline bool BoolMatHasSize(matlab::mex::ArgumentList& args, const int idx,
                             const int rows, const int cols) {
  return IsArgumentBoolMat(args, idx) and
         (args[idx].getDimensions().size() == 2) and
         (args[idx].getDimensions()[0] == rows) and
         (args[idx].getDimensions()[1] == cols);
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
void ValidateArgumentIsBoolMatWithSize(matlab::mex::ArgumentList& args,
                                         const int idx, const int rows,
                                         const int cols) {
  ValidateArgumentSizeOfAtLeast(args, idx);
  ValidateArgumentIsBoolMat(args, idx);
  if (not BoolMatHasSize(args, idx, rows, cols)) {
    if (args[idx].getDimensions().size() != 2) {
      throw std::runtime_error(
          "Expected 2 dimensions, but got " +
          std::to_string(args[idx].getDimensions().size()));
    }
    throw std::runtime_error(
        "Expected [" + std::to_string(rows) + " x " + std::to_string(cols) +
        "] bool matrix but got [" +
        std::to_string(args[idx].getDimensions()[0]) + " x " +
        std::to_string(args[idx].getDimensions()[1]) + "] instead");
  }
}
inline double ValidateAndLoadArgumentSingleDoubleValue(
    matlab::mex::ArgumentList& args, const int idx) {
  ValidateArgumentIsDoubleMatWithSize(args, idx, 1, 1);
  return double{args[idx][0][0]};
}
inline bool ValidateAndLoadArgumentSingleBoolValue(
    matlab::mex::ArgumentList& args, const int idx) {
  ValidateArgumentIsBoolMatWithSize(args, idx, 1, 1);
  return bool{args[idx][0][0]};
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
  constexpr int kNumExpectedInputs = 11;
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
  ret.rover_state_ = LoadRoverStateFromMatlab(inputs[kStateVarInputNum]);
  ret.waypoint_global_no_mirror_ =
      LoadWaypointFromMatlab(inputs[kWaypointInputNum]);
  ret.always_risky_obs_ =
      ProcessAlwaysRiskyInput(inputs[kAlwaysRiskyObsInputNum]);
  ret.always_non_risky_obs_ =
      ProcessDynObs(inputs[kAlwaysNonRiskyObsInputNum],
                    roahm::DynObs::ObstacleType::kStaticBoundary);
  ret.maybe_risky_obs_ = ProcessChooseable(inputs[kChooseableObsInputNum]);
  return ret;
}

class MexFunction : public matlab::mex::Function {
private:
  std::string last_frs_fpath_;
  ::roahm::RiskRtd risk_rtd_;

public:
  MexFunction() : last_frs_fpath_{"/data/cpp_processed_CUDA_Highway_frs.frs.bin"}, risk_rtd_{last_frs_fpath_} {}
  void operator()(matlab::mex::ArgumentList outputs,
                  matlab::mex::ArgumentList inputs) {
    using ::roahm::GetDeltaS;
    using ::roahm::Tick;
    try {
      matlab::data::ArrayFactory factory;
      matlab::data::TypedArray<double> arr_out =
          factory.createArray<double>({8});

      auto risk_inputs{ValidateAndLoadInputs(outputs, inputs)};
      const std::string new_frs_fpath{ValidateAndLoadSingleString(inputs, 5)};
      const double risk_threshold{ValidateAndLoadArgumentSingleDoubleValue(inputs, 6)};
      const bool check_mirrors{ValidateAndLoadArgumentSingleBoolValue(inputs, 7)};
      const bool use_waypoint_modification_heuristic{ValidateAndLoadArgumentSingleBoolValue(inputs, 8)};
      const bool use_selection_heuristic{ValidateAndLoadArgumentSingleBoolValue(inputs, 9)};
      const bool use_left_cost_function{ValidateAndLoadArgumentSingleBoolValue(inputs, 10)};
      //fmt::print("[LOADED RISK THRESHOLD] {}\n", risk_threshold);
      if (new_frs_fpath != last_frs_fpath_) {
        //fmt::print("[SETTING NEW FRS FPATH] {}\n", new_frs_fpath);
	last_frs_fpath_ = new_frs_fpath;
        //fmt::print("[LOADING FRS @ FPATH] {}\n", new_frs_fpath);
	risk_rtd_ = ::roahm::RiskRtd{last_frs_fpath_};
        //fmt::print("[LOADED FRS @ FPATH] {}\n", new_frs_fpath);
      }
      // fmt::print("[MEX] Ran Planning Iteration\n");
      // fmt::print("[MEX] RISK RTD Took: {}\n",
      // risk_outputs.total_time_seconds_);

      //fmt::print("[RRM] Loading FRS\n");
      const auto frs_full{roahm::LoadFrsBinary(last_frs_fpath_)};
      //fmt::print("[RRM] Loaded FRS\n");
      // fmt::print("Num Non Risky:   {}\n",
      //            risk_inputs.always_non_risky_obs_.GetNumObs());
      // fmt::print("Num Risky:       {}\n", risk_inputs.always_risky_obs_.size());
      // fmt::print("Num Maybe Risky: {}\n", risk_inputs.maybe_risky_obs_.size());

      // fmt::print("[RMM] Pre Planning ITER\n");
      const auto t0 = Tick();
      std::vector<roahm::RiskRtdPlanningOutputsIpopt> risk_outputs_original;
      risk_outputs_original = risk_rtd_.RunPlanningIteration<::roahm::risk_rtd_ipopt_problem::RiskRtdIpoptProblem>(risk_inputs, true, risk_threshold, check_mirrors, use_waypoint_modification_heuristic, use_left_cost_function);
      auto risk_outputs = risk_outputs_original;
      const auto t1 = Tick();

      const double delta_y_wp{risk_inputs.waypoint_global_no_mirror_.y_ -
                              risk_inputs.rover_state_.y_};
      const bool is_waypoint_in_same_lane{std::abs(delta_y_wp) <= 0.1};
      {
        int color_idx{0};
        int idx{0};
        int success_count{0};
        // fmt::print("[RRM] Iterating\n");
        for (const auto& full_opt_inf : risk_outputs) {
          ++idx;
          success_count += full_opt_inf.found_feasible_;
          const auto& curr_frs = frs_full.GetVehrs(full_opt_inf.sel_inf_);
          if (true or (not full_opt_inf.found_feasible_)) {
//            fmt::print(R"""([ RES ] Problem {} / {} (1-idx)
//  Mirror:        {}
//  ManuType:      {}
//  Feasible:      {}
//  Solve Success: {}  
//  Status:        {}
//  Param:         {} in [{}, {}]
//  Cost:          {}
//  Solver Status: {}
//)""",
//                       idx, risk_outputs.size(), full_opt_inf.sel_inf_.mirror_,
//                       ToString(full_opt_inf.sel_inf_.manu_type_),
//                       full_opt_inf.found_feasible_,
//                       full_opt_inf.solve_succeeded_,
//                       ::roahm::ToString(full_opt_inf.ipopt_status_),
//                       ::roahm::StdToStringOpt(full_opt_inf.final_param_),
//                       curr_frs.GetTrajParamMin(), curr_frs.GetTrajParamMax(),
//                       ::roahm::StdToStringOpt(full_opt_inf.final_cost_),
//                       ::roahm::ToStringOpt<Ipopt::SolverReturn>(
//                           full_opt_inf.final_solver_status_));

            const double param_val_unmirrored =
                full_opt_inf.final_param_.value_or(curr_frs.GetTrajParamMin());
          }
        }
        fmt::print("{} successful / {}\n", success_count, risk_outputs.size());
        // fmt::print("waypoint_in_same_lane: {}\n", is_waypoint_in_same_lane);
        // fmt::print("global waypoint: {} {} {}\n",
        //            risk_inputs.waypoint_global_no_mirror_.x_,
        //            risk_inputs.waypoint_global_no_mirror_.y_,
        //            risk_inputs.waypoint_global_no_mirror_.h_);
        const auto local_tmp =
            risk_inputs.waypoint_global_no_mirror_.RelativeToEgoFrameNoMirror(
                risk_inputs.rover_state_.GetXYH());

        if (success_count > 0) {
          // Found feasible
          const double u0{risk_inputs.rover_state_.u_};
          const double h0{risk_inputs.rover_state_.GetHeading()};
          const auto rtd_out_comparison_operator =
              [u0, h0, delta_y_wp, is_waypoint_in_same_lane, use_selection_heuristic](
                  const ::roahm::RiskRtdPlanningOutputs& p0,
                  const ::roahm::RiskRtdPlanningOutputs& p1) -> bool {
            // p0 < p1  ===>  p0 has higher priority
            constexpr bool kPrioritizeP0{true};
            constexpr bool kPrioritizeP1{false};
            static_assert(kPrioritizeP0 != kPrioritizeP1);

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

            if (not use_selection_heuristic) {
              if (p0_cost <= p1_cost) {
                return kPrioritizeP0;
              }
              return kPrioritizeP1;
            }

            const ::roahm::PointXYH p0_location_local{
                p0.final_location_.value()};
            const ::roahm::PointXYH p1_location_local{
                p1.final_location_.value()};
            const double p0_mirror_mult{p0.sel_inf_.mirror_ ? -1.0 : 1.0};
            const double p1_mirror_mult{p1.sel_inf_.mirror_ ? -1.0 : 1.0};
            const double p0_delta_y{p0_location_local.y_ * p0_mirror_mult};
            const double p1_delta_y{p1_location_local.y_ * p1_mirror_mult};
            const double p0_signed_dist_to_wp_y{delta_y_wp - p0_delta_y};
            const double p1_signed_dist_to_wp_y{delta_y_wp - p1_delta_y};
            const double p0_abs_dist_to_wp_y{std::abs(p0_signed_dist_to_wp_y)};
            const double p1_abs_dist_to_wp_y{std::abs(p1_signed_dist_to_wp_y)};

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

	  // Enable for debugging
          // for (int i = 0; i < risk_outputs.size(); ++i) {
          //   const auto& p0 = risk_outputs.at(i);
          //   fmt::print(
          //       "[{}]\n Type: {}\n P: {} of [{}, {}]\n WPLMY: {}\n FLY: {}\n",
          //       i, ToString(p0.sel_inf_.manu_type_), p0.final_param_.value(),
          //       frs_full.GetVehrs(p0.sel_inf_).GetTrajParamMin(),
          //       frs_full.GetVehrs(p0.sel_inf_).GetTrajParamMax(),
          //       p0.waypoint_local_frame_.y_, p0.final_location_.value().y_);
          // }
          if (not risk_outputs.front().found_feasible_) {
            throw std::runtime_error(
                "Sorted, success count > 0, but front doesn't have feasible?");
          }

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
          constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
          arr_out[0] = true;
          arr_out[1] = min_param_no_mirror.value_or(kNaN);
          arr_out[2] = min_mirror;
          arr_out[3] = min_manu_type_int;
          arr_out[4] = min_u0.value_or(kNaN);
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
          }
          if (::roahm::IsSpd(
                  min_manu_type.value_or(::roahm::ManuType::kNone))) {
            fmt::print("Successfully got outputs in {} sec, spd {}mps\n",
                       GetDeltaS(t1, t0),
                       min_param_no_mirror.value_or(kNaN));
          } else {
            fmt::print("Successfully got outputs in {} sec, spd {}mps\n",
                       GetDeltaS(t1, t0),
                       risk_inputs.rover_state_.u_);
          }
        } else {
          fmt::print("No successes\n");
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
