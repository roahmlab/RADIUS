#ifndef ROAHM_FRS_LOADER_HPP_
#define ROAHM_FRS_LOADER_HPP_

#include <fstream>
#include <string> // for string
#include <vector> // for vector

#include "common.hpp" // for IndexT
#include "frs_mega.hpp"
#include "frs_total.hpp"
#include "manu_type.hpp" // for ManuType
#include "rover_state.hpp"
#include "simple_util.hpp" // for Interval
#include "vehrs.hpp"

/// @file frs_loader.hpp Contains functions to load Forward Reachable Sets
/// (FRSes) from files

namespace roahm {

FrsTotal LoadFrsBinary(const std::string& file_name);

template <typename... Args> struct is_stl_vec : std::false_type {};
template <typename... Args>
struct is_stl_vec<std::vector<Args...>> : std::true_type {};

template <typename T> bool NearEq(const T& v0, const T& v1, double tol);

template <typename T>
bool VecsEq(const std::vector<T>& v0, const std::vector<T>& v1) {
  if (v0.size() != v1.size()) {
    std::cerr << "Vectors have different sizes, " << v0.size()
              << " != " << v1.size() << std::endl;
    return false;
  }
  for (std::size_t i = 0; i < v0.size(); ++i) {
    if constexpr (is_stl_vec<T>::value) {
      if (not VecsEq(v0.at(i), v1.at(i))) {
        return false;
      }
    } else {
      if (not(v0.at(i) == v1.at(i))) {
        if constexpr (std::is_same<T, double>::value ||
                      std::is_same<T, float>::value ||
                      std::is_same<T, int>::value ||
                      std::is_same<T, IndexT>::value ||
                      std::is_same<T, std::size_t>::value) {
          std::cout << "Vectors differ at index " << i << ", " << v0.at(i)
                    << " != " << v1.at(i) << std::endl;
        } else {
          std::cout << "Vectors differ at index " << i << std::endl;
        }
        return false;
      }
    }
  }
  return true;
}

template <typename T>
bool VecsEq(const std::vector<T>& v0, const std::vector<T>& v1,
            std::string error_info) {
  if (not VecsEq(v0, v1)) {
    std::cerr << "VecsEq Returned False: " << error_info << std::endl;
    return false;
  }
  return true;
}

template <typename T>
bool VecsNearEq(const std::vector<T>& v0, const std::vector<T>& v1,
                double tol) {
  if (v0.size() != v1.size()) {
    // TODO add this back
    // std::cerr << "Vectors have different sizes" << std::endl;
    std::cerr << "Vectors have different sizes, " << v0.size()
              << " != " << v1.size() << std::endl;
    return false;
  }
  for (std::size_t i = 0; i < v0.size(); ++i) {
    if constexpr (is_stl_vec<T>::value) {
      if (not VecsNearEq(v0.at(i), v1.at(i), tol)) {
        return false;
      }
    } else {
      if (not NearEq<T>(v0.at(i), v1.at(i), tol)) {
        // TODO clean this up
        if constexpr (std::is_same<T, double>::value ||
                      std::is_same<T, float>::value ||
                      std::is_same<T, int>::value ||
                      std::is_same<T, IndexT>::value ||
                      std::is_same<T, std::size_t>::value) {
          std::cout << "Vectors differ at index " << i << ", " << v0.at(i)
                    << " != " << v1.at(i) << std::endl;
          std::cout << "Delta: " << v0.at(i) - v1.at(i) << std::endl;
          std::cout << "Tol: " << tol << std::endl;
        } else {
          std::cout << "Vectors differ at index " << i << std::endl;
        }
        return false;
      }
    }
  }
  return true;
}

template <typename T>
bool VecsNearEq(const std::vector<T>& v0, const std::vector<T>& v1, double tol,
                std::string error_info) {
  if (not VecsNearEq(v0, v1, tol)) {
    std::cerr << "VecsNearEq Returned False: " << error_info << std::endl;
    return false;
  }
  return true;
}

template <> bool NearEq(const double& v0, const double& v1, double tol);
template <> bool NearEq(const Sliceable& v0, const Sliceable& v1, double tol);
template <>
bool NearEq(const ZonoSliceInfo& v0, const ZonoSliceInfo& v1, double tol);

template <> bool NearEq(const Interval& v0, const Interval& v1, double tol);

template <> bool NearEq(const Vehrs& v0, const Vehrs& v1, double tol);

template <> bool NearEq(const FrsMega& v0, const FrsMega& v1, double tol);

} // namespace roahm
#endif // ROAHM_FRS_LOADER_HPP_
