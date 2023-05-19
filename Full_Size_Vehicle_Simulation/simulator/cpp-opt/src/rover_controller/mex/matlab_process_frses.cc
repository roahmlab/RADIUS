// Load and pre-process FRSes from MAT file and store them in C++ data
// structures. This is intended to be run offline.
#include <algorithm>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "MatlabDataArray/ArrayType.hpp"
#include "MatlabDataArray/String.hpp"
#include "MatlabDataArray/TypedArray.hpp"
#include "frs_io.hpp"
#include "frs_loader.hpp"
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "process_utils.hpp"

std::string MatlabTypeToStr(matlab::data::ArrayType type);
::roahm::CudaInfo
GetCudaInfoFromMat(const matlab::data::StructArray& cuda_frs_mat_info);

int GetIntFromSingleEltMat(const Eigen::MatrixXd& matrix) {
  if (matrix.rows() != 1 || matrix.cols() != 1) {
    throw std::runtime_error("Matrix is not a single element.");
  }
  return matrix(0, 0);
}

std::vector<double> GetDoubleVectorFromRowMat(const Eigen::MatrixXd& matrix) {
  if (matrix.rows() != 1) {
    throw std::runtime_error("Matrix must have 1 row.");
  }
  if (matrix.cols() == 0) {
    throw std::runtime_error("Matrix must have at least 1 column.");
  }
  std::vector<double> vec;
  vec.reserve(matrix.cols());
  for (int i = 0; i < matrix.cols(); ++i) {
    vec.push_back(matrix(0, i));
  }
  return vec;
}
std::vector<int> GetIntVectorFromRowMat(const Eigen::MatrixXd& matrix) {
  if (matrix.rows() != 1) {
    throw std::runtime_error("Matrix must have 1 row.");
  }
  if (matrix.cols() == 0) {
    throw std::runtime_error("Matrix must have at least 1 column.");
  }
  std::vector<int> vec;
  vec.reserve(matrix.cols());
  for (int i = 0; i < matrix.cols(); ++i) {
    vec.push_back(matrix(0, i));
  }
  return vec;
}

template <typename InputType> auto GenLambda(InputType& bind_var);

template <> auto GenLambda(int& bind_var) {
  return [&bind_var](const Eigen::MatrixXd& mat) -> void {
    bind_var = GetIntFromSingleEltMat(mat);
  };
}
template <> auto GenLambda(std::vector<double>& bind_var) {
  return [&bind_var](const Eigen::MatrixXd& mat) -> void {
    bind_var = GetDoubleVectorFromRowMat(mat);
  };
}
template <> auto GenLambda(std::vector<int>& bind_var) {
  return [&bind_var](const Eigen::MatrixXd& mat) -> void {
    bind_var = GetIntVectorFromRowMat(mat);
  };
}

struct ExpectedField {
  const std::string name_;
  std::optional<Eigen::MatrixXd> in_value_;
  std::function<void(const Eigen::MatrixXd&)> converter_;

  template <typename T>
  ExpectedField(const std::string& name, T& bind_var)
      : name_{name}, in_value_{std::nullopt}, converter_{GenLambda(bind_var)} {}
};
::roahm::CudaInfo
GetCudaInfoFromMat(const matlab::data::StructArray& cuda_frs_mat_info) {
  const std::string grid_size_str = "grid_size";
  const std::string num_zono_str = "num_zono";
  const std::string grid_x0_str = "grid_x0";
  const std::string grid_y0_str = "grid_y0";
  const std::string grid_dx_str = "grid_dx";
  const std::string grid_dy_str = "grid_dy";
  const std::string g_u0_x_str = "g_u0_x";
  const std::string g_u0_y_str = "g_u0_y";
  const std::string g_v0_x_str = "g_v0_x";
  const std::string g_v0_y_str = "g_v0_y";
  const std::string g_r0_x_str = "g_r0_x";
  const std::string g_r0_y_str = "g_r0_y";
  const std::string block_inzono_list_str = "block_inzono_list";
  const std::string rot_angle_str = "rot_angle";
  const std::string g_p_x_str = "g_p_x";
  const std::string cg_p_str = "cg_p";
  const std::string manu_type_str = "manu_type";
  const std::string cg_u0v0r0_str = "cg_u0v0r0";
  const std::string t_range_str = "t_range";

  ::roahm::CudaInfo cuda_info;
  std::vector<ExpectedField> expected_fields{
      {grid_size_str, cuda_info.grid_size_},
      {num_zono_str, cuda_info.num_zono_},
      {grid_x0_str, cuda_info.grid_x0_},
      {grid_y0_str, cuda_info.grid_y0_},
      {grid_dx_str, cuda_info.grid_dx_},
      {grid_dy_str, cuda_info.grid_dy_},
      {g_u0_x_str, cuda_info.g_u0_x_},
      {g_u0_y_str, cuda_info.g_u0_y_},
      {g_v0_x_str, cuda_info.g_v0_x_},
      {g_v0_y_str, cuda_info.g_v0_y_},
      {g_r0_x_str, cuda_info.g_r0_x_},
      {g_r0_y_str, cuda_info.g_r0_y_},
      {block_inzono_list_str, cuda_info.block_inzono_list_},
      {rot_angle_str, cuda_info.rot_angle_},
      {g_p_x_str, cuda_info.g_p_x_},
      {cg_p_str, cuda_info.cg_p_},
      {manu_type_str, cuda_info.manu_type_},
      {cg_u0v0r0_str, cuda_info.cg_u0v0r0_},
      {t_range_str, cuda_info.t_range_}};

  const int num_fields = cuda_frs_mat_info.getNumberOfFields();
  const auto num_expected_fields = expected_fields.size();
  if (num_fields != num_expected_fields) {
    throw std::runtime_error(
        "Expected CUDA struct to have " + std::to_string(num_expected_fields) +
        " fields, but received " + std::to_string(num_fields));
  }
  for (const auto& field_name : cuda_frs_mat_info.getFieldNames()) {
    const std::string field_name_str = std::string(field_name);
    for (auto& expected_field : expected_fields) {
      if (field_name_str == expected_field.name_) {
        const auto& field_val = cuda_frs_mat_info[0][field_name];
        const auto& field_type = field_val.getType();
        if (field_type != matlab::data::ArrayType::DOUBLE) {
          throw std::runtime_error(
              "Expected field " + expected_field.name_ + " to be of type " +
              MatlabTypeToStr(matlab::data::ArrayType::DOUBLE) +
              ", but received " + MatlabTypeToStr(field_type));
        }
        const auto& dbl_arr = matlab::data::TypedArray<double>(field_val);
        const int num_dims = dbl_arr.getDimensions().size();
        constexpr int kExpectedNumDims = 2;
        if (num_dims != kExpectedNumDims) {
          throw std::runtime_error(
              "Expected field " + expected_field.name_ + " to have " +
              std::to_string(kExpectedNumDims) + " dimensions, but received " +
              std::to_string(num_dims));
        }
        const auto& field_type_str = MatlabTypeToStr(field_type);
        const int n_rows = dbl_arr.getDimensions()[0];
        const int n_cols = dbl_arr.getDimensions()[1];
        auto dbl_cp = dbl_arr;
        const auto uniq_ptr = dbl_cp.release();
        expected_field.in_value_ =
            Eigen::Map<Eigen::MatrixXd>(uniq_ptr.get(), n_rows, n_cols);
      }
    }
  }

  for (const auto& expected_field : expected_fields) {
    if (not expected_field.in_value_.has_value()) {
      throw std::runtime_error("Expected CUDA struct to have field " +
                               expected_field.name_);
    }

    expected_field.converter_(expected_field.in_value_.value());
  }

  return cuda_info;
}

/// @brief Convert a MATLAB array type to a string.
/// @param type the MATLAB array type to convert.
/// @return a string representation of the MATLAB array type.
std::string MatlabTypeToStr(matlab::data::ArrayType type) {
  switch (type) {
  case matlab::data::ArrayType::DOUBLE:
    return "double";
  case matlab::data::ArrayType::SINGLE:
    return "single";
  case matlab::data::ArrayType::INT8:
    return "int8";
  case matlab::data::ArrayType::UINT8:
    return "uint8";
  case matlab::data::ArrayType::INT16:
    return "int16";
  case matlab::data::ArrayType::UINT16:
    return "uint16";
  case matlab::data::ArrayType::INT32:
    return "int32";
  case matlab::data::ArrayType::UINT32:
    return "uint32";
  case matlab::data::ArrayType::INT64:
    return "int64";
  case matlab::data::ArrayType::UINT64:
    return "uint64";
  case matlab::data::ArrayType::CHAR:
    return "char";
  case matlab::data::ArrayType::LOGICAL:
    return "logical";
  case matlab::data::ArrayType::COMPLEX_DOUBLE:
    return "complex double";
  case matlab::data::ArrayType::COMPLEX_SINGLE:
    return "complex single";
  case matlab::data::ArrayType::COMPLEX_INT8:
    return "complex int8";
  case matlab::data::ArrayType::COMPLEX_UINT8:
    return "complex uint8";
  case matlab::data::ArrayType::COMPLEX_INT16:
    return "complex int16";
  case matlab::data::ArrayType::COMPLEX_UINT16:
    return "complex uint16";
  case matlab::data::ArrayType::COMPLEX_INT32:
    return "complex int32";
  case matlab::data::ArrayType::COMPLEX_UINT32:
    return "complex uint32";
  case matlab::data::ArrayType::COMPLEX_INT64:
    return "complex int64";
  case matlab::data::ArrayType::COMPLEX_UINT64:
    return "complex uint64";
  case matlab::data::ArrayType::CELL:
    return "cell";
  case matlab::data::ArrayType::STRUCT:
    return "struct";
  // case matlab::data::ArrayType::FUNCTION_HANDLE:
  //     return "function handle";
  case matlab::data::ArrayType::OBJECT:
    return "object";
  case matlab::data::ArrayType::UNKNOWN:
    return "unknown specified";
  default:
    return "unknown unspecified";
  }
  return "UKN";
}

/// @brief Find nonzero elements in an iterable container.
/// @tparam T the type of the container.
/// @param container the container to search.
/// @param offset the offset to add to the indices of the nonzero elements.
/// @return a vector of indices of nonzero elements, relative to the begin
/// iterator.
template <typename T> std::vector<int> FindNZ(const T& container, int offset) {
  std::vector<int> nnz;
  int ctr = offset;
  for (const auto& val : container) {
    if (val != 0) {
      nnz.push_back(ctr);
    }
    ++ctr;
  }
  return nnz;
}

template <typename T> std::vector<int> FindNZ(const T& container) {
  return FindNZ(container, 0);
}

/// @brief Checks whether every element in a vector is unique.
/// @param vec_in the vector whose elements are being checked.
/// @return True if all elements are unique, false otherwise.
bool OnlyUniqueElements(const std::vector<int>& vec_in) {
  auto vec = vec_in;
  std::sort(vec.begin(), vec.end());
  const auto last = std::unique(vec.begin(), vec.end());
  if (last != vec.end()) {
    return false;
  }
  return true;
}

/// @brief Swap the columns of a matrix.
/// @tparam T the type of the matrix.
/// @param mat the matrix whose columns are being swapped.
/// @param from the original indices of the columns to swap
/// @param to the new indices of the columns to swap, corresponding to the
/// original indices.
template <typename T>
void SwapCols(T& mat, std::vector<int> from, std::vector<int> to) {
  if (from.size() != to.size()) {
    throw std::runtime_error(
        "from and to vectors must be the same size to swap columns");
  }
  if (not OnlyUniqueElements(to)) {
    throw std::runtime_error("to vector must contain unique elements");
  }
  if (not OnlyUniqueElements(from)) {
    throw std::runtime_error("to vector must contain unique elements");
  }

  for (int i = 0; i < from.size(); ++i) {
    const auto to_idx = to.at(i);
    const auto from_idx = from.at(i);

    const Eigen::VectorXd to_val = mat.col(to_idx);
    const Eigen::VectorXd from_val = mat.col(from_idx);

    mat.col(to_idx) = from_val;
    mat.col(from_idx) = to_val;

    for (int j = i + 1; j < from.size(); ++j) {
      if (from.at(j) == to_idx) {
        from.at(j) = from_idx;
      }
    }
  }
}

std::string
LoadOutputFpathFromMatlabArg(const matlab::data::Array& out_fpath_matlab_str) {
  const auto out_fpath_matlab_str_type = out_fpath_matlab_str.getType();
  if (out_fpath_matlab_str_type != matlab::data::ArrayType::MATLAB_STRING) {
    throw std::runtime_error(
        "Expected output fpath specified in input to be a MATLAB string, "
        "but received other type.");
  }
  const auto str_vec = matlab::data::StringArray{out_fpath_matlab_str};
  if (str_vec.isEmpty()) {
    throw std::runtime_error("No output string provided");
  }
  const auto num_strs = std::distance(str_vec.begin(), str_vec.end());
  if (num_strs != 1) {
    throw std::runtime_error("Expected one string, not " +
                             std::to_string(num_strs));
  }
  const std::string out_fpath_str{str_vec[0]};
  std::cout << "out_fpath_str: " << out_fpath_str << std::endl;
  return out_fpath_str;
}

// using namespace matlab::data;
// using matlab::mex::ArgumentList;
//
template <typename T> class MexArray2d {
public:
  matlab::data::TypedArray<T> array_;
  MexArray2d(matlab::data::TypedArray<T>& array) : array_{std::move(array)} {}

  T& operator()(int i, int j) {
    return *(array_.begin() + (i + j * array_.getDimensions()[0]));
  }

  const T& operator()(int i, int j) const {
    return *(array_.begin() + (i + j * array_.getDimensions()[0]));
  }

  const int NumRows() const { return array_.getDimensions()[0]; }

  const int NumCols() const { return array_.getDimensions()[1]; }

  Eigen::MatrixXd ToEigenCopy() const {
    // Construct Eigen matrix from Matlab array
    const auto num_rows = array_.getDimensions()[0];
    const auto num_cols = array_.getDimensions()[1];
    Eigen::MatrixXd eigen_mat(num_rows, num_cols);
    // Fill matrix with Matlab array
    for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < num_cols; j++) {
        eigen_mat(i, j) = operator()(i, j);
      }
    }
    return eigen_mat;
  }
};
// Mex C++ to load a matrix passed by MATLAB
class MexFunction : public matlab::mex::Function {
public:
  auto GetTime() -> decltype(std::chrono::high_resolution_clock::now()) {
    return std::chrono::high_resolution_clock::now();
  }
  static double GetDeltaS(std::chrono::high_resolution_clock::time_point t2,
                          std::chrono::high_resolution_clock::time_point t1) {
    return static_cast<double>(
               std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)
                   .count()) /
           1.0e9;
  }

  template <typename T>
  void ProcessZono(const T& zono_elt, const bool is_spd, const bool is_lan,
                   ::roahm::Vehrs& vehrs,
                   std::optional<::roahm::Interval>& u0_interval) {
    if (zono_elt.getType() != matlab::data::ArrayType::DOUBLE) {
      throw std::runtime_error("Zono is of type " +
                               MatlabTypeToStr(zono_elt.getType()) +
                               " instead of DOUBLE");
    }
    const auto& zono_dbl = matlab::data::TypedArray<double>(zono_elt);
    const auto n_rows = zono_dbl.getDimensions()[0];
    const auto n_cols = zono_dbl.getDimensions()[1];

    if (n_rows != 20) {
      throw std::runtime_error("Zono has " + std::to_string(n_rows) +
                               " rows instead of 20");
    }

    if (n_cols < 2) {
      throw std::runtime_error("Zono has " + std::to_string(n_cols) +
                               " columns instead of at least 2");
    }

    auto cp = zono_dbl;
    auto uniq_ptr = cp.release();
    auto stride = Eigen::Stride<0, 0>();
    Eigen::Map<Eigen::MatrixXd> eigen_mat(uniq_ptr.get(), n_rows, n_cols,
                                          stride);
    std::vector<int> sliceable_param_dims;
    {
      // Remember, these are 0-indexed, not 1-indexed like in
      // MATLAB
      const int p_u_idx = 10;
      const int p_y_idx = 11;
      if (is_spd) {
        sliceable_param_dims.push_back(p_u_idx);
      } else {
        sliceable_param_dims.push_back(p_y_idx);
      }
    }

    std::vector<int> sliceable_param_idxs;
    for (const auto& param_idx : sliceable_param_dims) {
      const auto nz_param_vals =
          FindNZ(eigen_mat.row(param_idx).rightCols(n_cols - 1), 1);
      if (nz_param_vals.empty()) {
        throw std::runtime_error("No non-zero values for parameter " +
                                 std::to_string(param_idx));
      } else if (nz_param_vals.size() > 1) {
        throw std::runtime_error(
            "Zono has more than one nonzero parameter value on "
            "index " +
            std::to_string(param_idx));
      } else if (nz_param_vals.size() == 1) {
        sliceable_param_idxs.push_back(nz_param_vals.at(0));
      }
    }

    const int num_params_in_zono = sliceable_param_idxs.size();
    if (num_params_in_zono <= 0) {
      throw std::runtime_error("Zono has no sliceable parameters");
    }

    // TODO remove this if we are doing multiple parameters
    if (num_params_in_zono > 1) {
      throw std::runtime_error("Zono has more than one sliceable parameter");
    }

    std::vector<int> sliceable_param_to_idxs(num_params_in_zono);
    for (int i = 0; i < num_params_in_zono; i++) {
      sliceable_param_to_idxs.at(i) = n_cols - 1 - i;
    }
    std::reverse(sliceable_param_to_idxs.begin(),
                 sliceable_param_to_idxs.end());
    {
      // TODO cleanup
      const int slc_pre_col = sliceable_param_idxs.at(0);
      const int slc_row = sliceable_param_dims.at(0);
      const double pre_val = eigen_mat(slc_row, slc_pre_col);
      if (std::abs(pre_val) < 1.0e-9) {
        throw std::runtime_error("Zero Slice value: " +
                                 std::to_string(pre_val));
      }
      ::roahm::MoveColsToEnd(eigen_mat, sliceable_param_idxs);
    }

    constexpr int kUDim = 6;
    constexpr int kVDim = 7;
    constexpr int kRDim = 8;

    constexpr std::array<int, 3> kStateSliceableDimensions{kUDim, kVDim, kRDim};
    std::vector<int> state_sliceable_gen_idxs;
    for (const auto& dim : kStateSliceableDimensions) {
      // Minus one to account for the fact that the first column
      // is the center
      const auto nz_param_vals =
          FindNZ(eigen_mat.row(dim).rightCols(n_cols - 1), 1);
      if (nz_param_vals.size() > 1) {
        throw std::runtime_error(
            "Zono has more than one nonzero generator value on "
            "dimension " +
            std::to_string(dim));
      } else if (nz_param_vals.empty()) {
        throw std::runtime_error("Zono has no nonzero generator value on "
                                 "dimension " +
                                 std::to_string(dim));
      }
      if (dim == kUDim) {
        const double center_u = eigen_mat(kUDim, 0);
        const double gen_u = eigen_mat(kUDim, nz_param_vals.at(0));
        const double max_u0 = center_u + std::abs(gen_u);
        const double min_u0 = center_u - std::abs(gen_u);
        if (u0_interval) {
          const double prev_min_u0 = (*u0_interval).Min();
          const double prev_max_u0 = (*u0_interval).Max();
          if (min_u0 != prev_min_u0 || max_u0 != prev_max_u0) {
            throw std::runtime_error("u0 interval is not consistent. New: [" +
                                     std::to_string(min_u0) + ", " +
                                     std::to_string(max_u0) + "] Old: [" +
                                     std::to_string(prev_min_u0) + ", " +
                                     std::to_string(prev_max_u0) + "]");
          }
        } else {
          u0_interval = ::roahm::Interval{min_u0, max_u0};
        }
      }
      // Plus one to account for center not being searched
      state_sliceable_gen_idxs.push_back(nz_param_vals.at(0));
    }

    std::vector<int> state_sliceable_to_idxs;
    for (int i = 0; i < kStateSliceableDimensions.size(); ++i) {
      const int matrix_last_col_idx = n_cols - 1;
      const int dist_from_end_of_matrix = i + num_params_in_zono;
      const int des_idx = matrix_last_col_idx - dist_from_end_of_matrix;
      state_sliceable_to_idxs.push_back(des_idx);
    }
    std::reverse(state_sliceable_to_idxs.begin(),
                 state_sliceable_to_idxs.end());
    {
      auto inner_portion =
          eigen_mat.leftCols(eigen_mat.cols() - sliceable_param_idxs.size());
      // TODO REMOVE std::cout << "MAT BEFORE: \n" << eigen_mat << std::endl;
      ::roahm::MoveColsToEnd(inner_portion, state_sliceable_gen_idxs);
      // TODO REMOVE std::cout << "MAT AFTER: \n" << eigen_mat << std::endl;
    }
    // TODO REMOVE SwapCols(eigen_mat, state_sliceable_gen_idxs,
    // state_sliceable_to_idxs);

    for (int i = 0; i < kStateSliceableDimensions.size(); ++i) {
      const auto dim = kStateSliceableDimensions.at(i);
      // Minus one to account for the fact that the first column
      // is the center
      const auto nz_param_vals =
          FindNZ(eigen_mat.row(dim).rightCols(n_cols - 1), 1);
      if (nz_param_vals.size() > 1) {
        throw std::runtime_error("[After Swap] Zono has more than one nonzero "
                                 "generator value on dimension " +
                                 std::to_string(dim));
      } else if (nz_param_vals.empty()) {
        throw std::runtime_error(
            "[After Swap] Zono has no nonzero generator value "
            "on dimension " +
            std::to_string(dim));
      }

      const int state_sliceable_idx = nz_param_vals.at(0);
      const int expected_idx = state_sliceable_to_idxs.at(i);
      if (state_sliceable_idx != expected_idx) {
        throw std::runtime_error("Expected state sliceable index " +
                                 std::to_string(state_sliceable_idx) +
                                 " dim [" + std::to_string(dim) + "] to be " +
                                 std::to_string(expected_idx));
      }
    }

    // Find time sliceable generator
    // NOTE: for multiple time generators, we use sum of absolute values.
    const int kTimeDim = n_rows - 1;
    const auto nz_time_gens =
        FindNZ(eigen_mat.row(kTimeDim).rightCols(n_cols - 1), 1);
    if (nz_time_gens.empty()) {
      throw std::runtime_error("Zono has no nonzero generator value on time "
                               "dimension");
    }
    // if (nz_time_gens.size() > 1) {
    //   std::cout << "Eigen Mat: " << eigen_mat.row(kTimeDim) << std::endl;
    //   throw std::runtime_error("Zono has more than one nonzero generator
    //   value on "
    //                            "time dimension");
    // }
    const double time_center_val = eigen_mat(kTimeDim, 0);
    // TODO ACTUAL set this to be correct
    constexpr bool kUsingSum = false;
    constexpr bool kUsingFirst = true;
    constexpr bool kUsingMax = false;
    static_assert(kUsingSum xor kUsingFirst xor kUsingMax,
                  "Must use one of sum, first, or max");
    static_assert(kUsingSum or kUsingMax or kUsingFirst,
                  "Must use either sum or max or first");

    double total_abs_gen_sum =
        kUsingSum ? 0.0 : (-std::numeric_limits<double>::infinity());

    auto sorted_nz_time_gens = nz_time_gens;
    std::sort(sorted_nz_time_gens.begin(), sorted_nz_time_gens.end());

    for (const auto& nz_time_gen_col_idx : sorted_nz_time_gens) {
      // TODO check this
      const auto curr_gen_val =
          std::abs(eigen_mat(kTimeDim, nz_time_gen_col_idx));
      if constexpr (kUsingSum) {
        // std::cout << "Adding " << curr_gen_val << " to total_abs_gen_sum" <<
        // std::endl;
        total_abs_gen_sum += curr_gen_val;
      } else if (kUsingMax) {
        total_abs_gen_sum = std::max(total_abs_gen_sum, curr_gen_val);
      } else if (kUsingFirst) {
        total_abs_gen_sum = curr_gen_val;
        break;
      }
    }
    // std::cout << "Total abs gen sum: " << total_abs_gen_sum << std::endl;

    // TODO very very sketch
    // TODO ACTUAL add this back
    // if (nz_time_gens.size() > 1) {
    //   total_abs_gen_sum = 0.01;
    // }
    const double time_max = time_center_val + total_abs_gen_sum;
    const double time_min = time_center_val - total_abs_gen_sum;
    // FrsZono<1> zono;
    //  TODO

    // Rows/dimensions of the zonotope correspond to the listed parameters
    constexpr int kXDim = 0;
    constexpr int kYDim = 1;
    constexpr int kHDim = 2;

    // Dimensionality of the decision variable
    constexpr int kDecisionVarDim = 1;

    const int curr_zono_idx = vehrs.GetNumZonos();
    const double zono_x_center = eigen_mat(kXDim, 0);
    const double zono_y_center = eigen_mat(kYDim, 0);
    const double zono_h_center = eigen_mat(kHDim, 0);

    vehrs.xy_centers_.push_back(zono_x_center);
    vehrs.xy_centers_.push_back(zono_y_center);
    vehrs.h_centers_.push_back(zono_h_center);

    // minus one for center
    const int left_offset = 1;
    const int right_offset = sliceable_param_dims.size();
    // + kStateSliceableDimensions.size()
    const auto eigen_non_slice_gens =
        eigen_mat.rightCols(n_cols - left_offset)
            .leftCols(n_cols - right_offset - left_offset);
    const int num_non_slice_gens = n_cols - left_offset - right_offset;
    if (num_non_slice_gens != eigen_non_slice_gens.cols()) {
      std::cout << "num_non_slice_gens [" << num_non_slice_gens
                << "] != eigen_non_slice_gens.cols() ["
                << eigen_non_slice_gens.cols() << "]" << std::endl;
      throw std::runtime_error(
          "num_non_slice_gens != eigen_non_slice_gens.cols()");
    }
    {
      if (vehrs.zono_sizes_.empty()) {
        std::cout << "[DBG] Right offset: " << right_offset << std::endl;
        std::cout << "[DBG] Left offset: " << left_offset << std::endl;
        std::cout << "[DBG] N Cols: " << n_cols << std::endl;
        std::cout << "[DBG] Num non slice gens: " << num_non_slice_gens
                  << std::endl;
      }
      int zono_gen_count = 0;
      int col_idx = left_offset;
      for (; col_idx < n_cols - right_offset; ++col_idx) {
        const double x_val = eigen_mat(kXDim, col_idx);
        const double y_val = eigen_mat(kYDim, col_idx);
        constexpr double kGenEps = 1.0e-5;
        if (std::abs(x_val) < kGenEps and std::abs(y_val) < kGenEps) {
          continue;
        }
        vehrs.zono_xy_.push_back(x_val);
        vehrs.zono_xy_.push_back(y_val);
        vehrs.zono_h_.push_back(eigen_mat(kHDim, col_idx));
        ++zono_gen_count;
      }
      if (col_idx != num_non_slice_gens + 1) {
        throw std::runtime_error("col_idx != num_non_slice_gens");
      }
      vehrs.zono_sizes_.push_back(zono_gen_count);
    }
    const ::roahm::Interval time_interval{time_min, time_max};
    std::cout << "[TDB] Time Interval: " << time_interval.Width() << "["
              << time_interval.Min() << ", " << time_interval.Max() << "]"
              << std::endl;
    vehrs.zono_time_intervals_.push_back(time_interval);
    // TODO slc_infos_
    //  vehrs.slc_x_.push_back(0.0);
    //  vehrs.slc_y_.push_back(0.0);
    //  vehrs.slc_xy_set_.push_back(false);
    // TODO t_eval_idx_
    const double des_time = is_lan ? 6.0 : 3.0;
    {
      double min_dist_to_des_time = std::numeric_limits<double>::max();
      int min_idx = -1;
      for (int i = 0; i < vehrs.zono_time_intervals_.size(); ++i) {
        const auto& t_int = vehrs.zono_time_intervals_.at(i);
        // TODO
        // const double curr_dist = t_int.DistanceTo(des_time);
        // TODO which to choose?
        const double curr_dist = std::abs(t_int.Midpoint() - des_time);
        if (curr_dist < min_dist_to_des_time) {
          min_dist_to_des_time = curr_dist;
          min_idx = i;
        }
      }
      vehrs.t_eval_idx_ = min_idx;
    }

    vehrs.slc_infos_.emplace_back();

    for (int i = 0; i < kStateSliceableDimensions.size(); ++i) {
      // Zero indexed
      const int sliceable_param_row = kStateSliceableDimensions.at(i);
      const int sliceable_param_col = state_sliceable_to_idxs.at(i);
      const double slc_center = eigen_mat(sliceable_param_row, 0);
      const double x = eigen_mat(kXDim, sliceable_param_col);
      const double y = eigen_mat(kYDim, sliceable_param_col);
      const double h = eigen_mat(kHDim, sliceable_param_col);
      const double slc_gen =
          eigen_mat(sliceable_param_row, sliceable_param_col);
      static_assert(kDecisionVarDim == 1,
                    "TODO k_rng isn't defined for kDecisionVarDim != 1");
      vehrs.slc_infos_.back().slc_vals_.emplace_back(
          sliceable_param_row + 1, slc_center, slc_gen, x, y, h);
    }

    for (int i = 0; i < sliceable_param_to_idxs.size(); ++i) {
      // Zero indexed
      const int sliceable_param_row = sliceable_param_dims.at(i);
      const int sliceable_param_col = sliceable_param_to_idxs.at(i);
      const double p_center = eigen_mat(sliceable_param_row, 0);
      const double x = eigen_mat(kXDim, sliceable_param_col);
      const double y = eigen_mat(kYDim, sliceable_param_col);
      const double h = eigen_mat(kHDim, sliceable_param_col);
      const double p_gen = eigen_mat(sliceable_param_row, sliceable_param_col);
      static_assert(kDecisionVarDim == 1,
                    "TODO k_rng isn't defined for kDecisionVarDim != 1");
      if (curr_zono_idx == 0) {
        // std::cout << "[EMPTY] vehrs.slc_xy_set_ is empty" << std::endl;
        vehrs.k_rng_ = p_gen;
      } else if (vehrs.k_rng_ != p_gen) {
        throw std::runtime_error("vehrs.k_rng_ [" +
                                 std::to_string(vehrs.k_rng_) + "] != p_gen [" +
                                 std::to_string(p_gen) + "]");
      }
      vehrs.slc_infos_.back().slc_vals_.emplace_back(sliceable_param_row + 1,
                                                     p_center, p_gen, x, y, h);
    }

    /*
    const ::roahm::Interval<double> time_interval{time_min, time_max};
    static_assert(kXDim + 1 == kYDim,
                  "need to change references if this does not hold");
    const int num_gen_cols = n_cols - 1 - sliceable_param_idxs.size() -
                             state_sliceable_gen_idxs.size();
    assert(num_gen_cols > 1 &&
           "removed most/all generator columns, something is probably wrong");
    const ::roahm::Zono2d xy_zono{eigen_mat.block(kXDim, 0, 2, 1),
                                  eigen_mat.block(kXDim, 1, 2, num_gen_cols)};
    const Eigen::Matrix<double, 2, kDecisionVarDim>
    decision_var_slice_generators_xy const Eigen::Matrix<double, 3,
    kDecisionVarDim> decision_var_xyh_generators const double center_h =
    eigen_mat(kHDim, 0); const UvrSliceable u_sliceable const UvrSliceable
    v_sliceable const UvrSliceable r_sliceable

            FrsZono<1>
            zono(time_interval, xy_zono, decision_var_slice_generators_xy,
                 decision_var_xyh_generators, center_h, u_sliceable,
                 v_sliceable, r_sliceable)
    */
  }

  template <typename T>
  ::roahm::Vehrs ProcessVehRS(const T& vehrs_cell_arr, const bool is_spd,
                              const bool is_lan,
                              std::optional<::roahm::Interval>& u0_interval) {
    ::roahm::Vehrs vehrs;
    int zono_ctr = 0;
    for (const auto& zono_elt : vehrs_cell_arr) {
      ProcessZono(zono_elt, is_spd, is_lan, vehrs, u0_interval);
      // std::cout << "k_rng: " << vehrs.k_rng_ << std::endl;
    }
    return vehrs;
  }

  void operator()(matlab::mex::ArgumentList outputs,
                  matlab::mex::ArgumentList inputs) {
    {
      // Verify Inputs
      constexpr int kNumExpectedInputs = 2;
      constexpr int kNumExpectedOutputs = 0;
      const int num_inputs = inputs.size();
      const int num_outputs = outputs.size();
      if (num_inputs != kNumExpectedInputs) {
        throw std::invalid_argument(
            "Expected " + std::to_string(kNumExpectedInputs) +
            " inputs, but got " + std::to_string(num_inputs));
      }
      if (num_outputs != kNumExpectedOutputs) {
        throw std::invalid_argument(
            "Expected " + std::to_string(kNumExpectedOutputs) +
            " outputs, but got " + std::to_string(num_outputs));
      }
    }

    const auto& top_level_arr = inputs[0];
    {
      const auto top_level_arr_type = top_level_arr.getType();
      if (top_level_arr_type != matlab::data::ArrayType::STRUCT) {
        throw std::runtime_error(
            "Expected input to be a struct, but received other type.");
      }
    }
    const auto& top_level_struct_arr = matlab::data::StructArray(top_level_arr);

    const std::string out_fpath = LoadOutputFpathFromMatlabArg(inputs[1]);
    // const matlab::data::String out_fpath_matlab_str_arr =
    // out_fpath_matlab_str; const std::string out_fpath =
    // static_cast<std::string>(matlab::data::StringArray(out_fpath_matlab_str));

    {
      const int num_fields = top_level_struct_arr.getNumberOfFields();
      if (num_fields != 1) {
        throw std::runtime_error(
            "Expected top level struct to have 1 field, but received " +
            std::to_string(num_fields));
      }

      auto field_names = top_level_struct_arr.getFieldNames();
      const std::string field_name = std::string(*field_names.begin());
      const std::string expected_m_mega_field_name = "M_mega";
      if (field_name != expected_m_mega_field_name) {
        throw std::runtime_error(
            "Expected top level struct to have field named " +
            expected_m_mega_field_name + ", but received " + field_name);
      }
    }

    const auto& all_megas = *top_level_struct_arr.begin();
    for (const auto& mega_it : all_megas) {
      const auto t = mega_it.getType();
      if (t != matlab::data::ArrayType::CELL) {
	std::cout << "Expected Mega to be CELL" << std::endl;
        return;
      }
      std::cout << "Mega_it: " << MatlabTypeToStr(t) << std::endl;
    }
    const auto& mega_cell_arr = matlab::data::CellArray(*all_megas.begin());
    std::cout << "Mega Cell Arr: " << mega_cell_arr.getNumberOfElements()
              << std::endl;
    int mega_ctr = 0;
    ::roahm::FrsTotal frs_total;
    for (const auto& single_mega : mega_cell_arr) {
      std::cout << "MEGA " << mega_ctr << std::endl;
      auto t = single_mega.getType();
      std::cout << "Val: " << MatlabTypeToStr(t) << std::endl;
      if (not(t == matlab::data::ArrayType::STRUCT)) {
        throw std::runtime_error("Expected MEGA to have type struct, but got " +
                                 MatlabTypeToStr(t));
      }

      const auto& single_mega_struct = matlab::data::StructArray(single_mega);
      auto field_names = single_mega_struct.getFieldNames();
      auto& mega_curr = frs_total.megas_.emplace_back();
      auto& au = mega_curr.au_;
      auto& dir = mega_curr.dir_;
      auto& lan = mega_curr.lan_;
      std::optional<::roahm::Interval> u0_interval;
      for (const auto field_name : field_names) {
        const auto& field_val = single_mega_struct[0][field_name];
        const auto field_val_type = field_val.getType();
        const std::string field_name_str = std::string(field_name);
        std::cout << "Field: " << field_name_str << std::endl;
        std::cout << " Type: " << MatlabTypeToStr(field_val_type) << std::endl;
        const bool is_spd = field_name_str == "Au";
        const bool is_lan = field_name_str == "lan";
        const bool is_dir = field_name_str == "dir";
        if (is_spd or is_lan or is_dir) {
          std::cout << "[DBG] FRS SET" << std::endl;
          if (field_val_type == matlab::data::ArrayType::CELL) {
            int num_cells = field_val.getNumberOfElements();
            std::cout << "Num Elts: " << num_cells << std::endl;
          } else {
            throw std::runtime_error("Expected field " + field_name_str +
                                     " to be a cell array, but received " +
                                     MatlabTypeToStr(field_val_type));
          }

          int frs_ctr = 0;
          const matlab::data::TypedArray<matlab::data::Array>&
              frs_individ_tarr = field_val;
          const auto& frs_manu_set_cell_arr =
              matlab::data::CellArray(frs_individ_tarr);
          std::cout << "frs_manu_set_cell_arr: "
                    << frs_manu_set_cell_arr.getNumberOfElements() << std::endl;
          for (auto manu_set_dim : frs_manu_set_cell_arr.getDimensions()) {
            std::cout << "manu_set_dim: " << manu_set_dim << std::endl;
          }
          for (const auto& frs_individ : frs_manu_set_cell_arr) {
            const auto frs_individ_type = frs_individ.getType();
            const auto frs_individ_type_str = MatlabTypeToStr(frs_individ_type);
            std::cout << "FRS [" << frs_ctr << "] = " << frs_individ_type_str
                      << std::endl;
            // TODO make this reasoning cleaner

            // Check if it's a double
            if (frs_individ_type == matlab::data::ArrayType::DOUBLE) {
              // TODO make sure it's either empty or full and not something
              // weird going on
              continue;
            }
            if (frs_individ_type != matlab::data::ArrayType::STRUCT) {
              throw std::runtime_error("Expected FRS structure, but received " +
                                       frs_individ_type_str);
            }

            const auto& frs_individ_struct =
                matlab::data::StructArray(frs_individ);

            {
              const std::string vehrs_save_field_str = "vehRS_save";
              const std::string cuda_info_field_str = "cuda_FRS";
              auto field_names = frs_individ_struct.getFieldNames();
              {
                struct ExpectedField {
                  std::string field_str_;
                  bool has_field_;
                };
                std::vector<ExpectedField> expected_fields{
                    {vehrs_save_field_str, false},
                    {cuda_info_field_str, false}};

                for (const auto& field_name : field_names) {
                  const auto field_name_str = std::string(field_name);
                  for (auto& expected_field : expected_fields) {
                    if (field_name_str == expected_field.field_str_) {
                      expected_field.has_field_ = true;
                    }
                  }
                }
                for (const auto& expected_field : expected_fields) {
                  if (!expected_field.has_field_) {
                    throw std::runtime_error(
                        "Expected field " + expected_field.field_str_ +
                        " in FRS structure, but did not find it");
                  }
                }
              }

              const matlab::data::MATLABFieldIdentifier vehrs_save_field_id{
                  vehrs_save_field_str};
              const matlab::data::MATLABFieldIdentifier cuda_frs_field_id{
                  cuda_info_field_str};
              ::roahm::CudaInfo cuda_info;
              {
                const auto& field_val =
                    frs_individ_struct[0][cuda_frs_field_id];
                const auto& field_type = field_val.getType();
                const auto& field_type_str = MatlabTypeToStr(field_type);
                std::cout << "[DBG] CUDA FRS FIELD" << std::endl;
                if (field_type != matlab::data::ArrayType::STRUCT) {
                  throw std::runtime_error(
                      "CheckFrs failed, cuda_FRS is not a cell");
                }
                const auto& cuda_frs_mat_struct =
                    matlab::data::StructArray(field_val);
                std::cout << "Pre CUDA INFO" << std::endl;
                cuda_info = GetCudaInfoFromMat(cuda_frs_mat_struct);
                std::cout << "Post CUDA INFO" << std::endl;
              }

              bool has_vehrs = false;
              {
                const auto& field_val =
                    frs_individ_struct[0][vehrs_save_field_id];
                const auto& field_type = field_val.getType();
                const auto& field_type_str = MatlabTypeToStr(field_type);
                if (field_type != matlab::data::ArrayType::CELL) {
                  throw std::runtime_error(
                      "CheckFrs failed, vehRS is not a cell");
                }
                has_vehrs = true;
                const matlab::data::TypedArray<matlab::data::Array>&
                    vehrs_typed_arr = field_val;
                const auto& vehrs_cell_arr =
                    matlab::data::CellArray(vehrs_typed_arr);
                auto vehrs_out{
                    ProcessVehRS(vehrs_cell_arr, is_spd, is_lan, u0_interval)};
                vehrs_out.cuda_info_ = cuda_info;
                if (is_spd) {
                  au.push_back(vehrs_out);
                } else if (is_dir) {
                  dir.push_back(std::vector<::roahm::Vehrs>{vehrs_out});
                } else if (is_lan) {
                  const int num_rows = frs_manu_set_cell_arr.getDimensions()[0];
                  const int num_cols = frs_manu_set_cell_arr.getDimensions()[1];
                  const int num_t_idxs =
                      frs_manu_set_cell_arr.getDimensions()[2];
                  if (lan.size() >= num_rows) {
                    throw std::runtime_error("lan.size() >= num_rows");
                  }
                  // TODO figure out what this should actually be
                  // if (lan.empty()) {
                  //}
                  lan.emplace_back(std::vector<::roahm::Vehrs>{vehrs_out});
                  // if ((lan.back().size() >= (num_cols * num_t_idxs))) {
                  //   lan.back.emplace_back();
                  // }
                  //  Ideally, this is how it should be?
                  //  if (lan.size() <= row) {
                  //    lan.emplace_back();
                  //  }
                  //  if (lan.back().size() <= col) {
                  //    lan.back().emplace_back();
                  //  }
                  //  if (lan.back().back().size() <= t_idx) {
                  //    lan.back().back().push_back(vehrs_out);
                  //  } else {
                  //    throw std::runtime_error("lan already has entry at [row,
                  //    col, t_idx]: [" +
                  //                             std::to_string(row) + ", " +
                  //                             std::to_string(col) + ", " +
                  //                             std::to_string(t_idx) + "]");
                  //  }
                }
              }

              if (not has_vehrs) {
                throw std::runtime_error(
                    "CheckFrs failed, vehRS_save not found");
              }
            }
            ++frs_ctr;
          }
        }
      }
      if (not u0_interval) {
        throw std::runtime_error("CheckFrs failed, u0_interval not found");
      } else {
        frs_total.u0_intervals_.push_back(*u0_interval);
      }
      ++mega_ctr;
    }

    if (not frs_total.u0_intervals_.empty()) {
      double min_u0_total = std::numeric_limits<double>::max();
      double max_u0_total = std::numeric_limits<double>::min();
      for (const auto& u0_interval : frs_total.u0_intervals_) {
        min_u0_total = std::min(min_u0_total, u0_interval.Min());
        max_u0_total = std::max(max_u0_total, u0_interval.Max());
      }
      frs_total.u0_min_ = min_u0_total;
      frs_total.u0_max_ = max_u0_total;
    } else {
      throw std::runtime_error("CheckFrs failed, u0_intervals_ is empty");
    }
    // TODO properly set this
    frs_total.alternating_au_ = false;
    frs_total.successful_ = true;

    if constexpr (true) {
      const std::string file_name{out_fpath};
      std::ofstream file_out(file_name, std::ios::binary);
      ::roahm::WriteToBinFile(file_out, frs_total);
      file_out.close();
    }

    return;
  }
};
