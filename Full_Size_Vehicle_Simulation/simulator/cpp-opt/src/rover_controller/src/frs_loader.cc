#include "frs_loader.hpp"

#include <algorithm> // for max, min_element
#include <cassert>   // for assert
#include <cstdlib>   // for abs
#include <fstream>   // for ifstream, ofstream
#include <iostream>  // for operator<<, ostream, cerr, istream, basic_...
#include <iterator>  // for distance
#include <tuple>
#include <type_traits> // for false_type, true_type

#include "frs_io.hpp"
#include "frs_loader.hpp"

namespace {
using IndexT = ::roahm::IndexT;
using Vehrs = ::roahm::Vehrs;
using FrsMega = ::roahm::FrsMega;
using Interval = ::roahm::Interval;

} // namespace

namespace roahm {

FrsTotal LoadFrsBinary(const std::string& file_name) {
  std::ifstream file_in(file_name, std::ios::binary);
  if (not file_in) {
    throw std::runtime_error("Could not open file at " + file_name);
  }
  const FrsTotal frs{ReadFromBinFile<FrsTotal>(file_in)};
  file_in.close();
  return frs;
}

SlicedInfo ZonoSliceInfo::Slice(double u, double v, double r) const {
  SlicedInfo ret;
  ret.lambda_ok_ = true;
  ret.slc_xy_set_ = false;
  for (const auto& sliceable : slc_vals_) {
    const auto sl_dim = sliceable.dim_;
    const auto sl_val = sliceable.slc_val_;
    // TODO what are these dims?
    if (sl_dim == 11 || sl_dim == 12) {
      ret.slc_x_ = sliceable.x_;
      ret.slc_y_ = sliceable.y_;
      ret.slc_h_ = sliceable.h_;
      ret.slc_val_ = sl_val;
      ret.center_slc_val_ = sliceable.center_val_;
      ret.slc_xy_set_ = true;
    } else {
      double slice_pt = 0.0;
      if (sl_dim == 7) {
        slice_pt = u;
      } else if (sl_dim == 8) {
        slice_pt = v;
      } else if (sl_dim == 9) {
        slice_pt = r;
      } else {
        assert(false);
      }
      constexpr double kLambdaEps = 1.0e-3;
      double slice_lambda = (slice_pt - sliceable.center_val_) / sl_val;
      double abs_lambda = std::abs(slice_lambda);
      if (abs_lambda > 1.0 && abs_lambda <= 1.0 + kLambdaEps) {
        if (slice_lambda < 0) {
          slice_lambda = -1.0;
        } else {
          slice_lambda = 1.0;
        }
      }
      abs_lambda = std::abs(slice_lambda);
      if (abs_lambda > 1.0) {
        std::cout << "Slice Point Out of bounds\n";
        std::cout << "slice_pt:              " << slice_pt << std::endl;
        std::cout << "sliceable.center_val_: " << sliceable.center_val_
                  << std::endl;
        std::cout << "sl_dim: " << sl_dim << std::endl;
        std::cout << "sl_val:    " << sl_val << std::endl;
        std::cout << "sl_lambda: " << slice_lambda << std::endl;
      }
      ret.lambda_ok_ &= (std::abs(slice_lambda) <= 1.0);
      ret.x_sliced_sum_ += slice_lambda * sliceable.x_;
      ret.y_sliced_sum_ += slice_lambda * sliceable.y_;
      ret.h_sliced_sum_ += slice_lambda * sliceable.h_;
    }
  }
  return ret;
}

bool ZonoSliceInfo::operator==(const ZonoSliceInfo& oth) const {
  return VecsEq(slc_vals_, oth.slc_vals_);
}

bool Sliceable::NearEqual(const Sliceable& other, double tol) const {
  // std::cout << "Sliceable::NearEqual" << std::endl;
  // std::cout << "[This] " << *this << std::endl;
  // std::cout << "[Oth] " << other << std::endl;
  return dim_ == other.dim_ && NearEq(center_val_, other.center_val_, tol) &&
         NearEq(slc_val_, other.slc_val_, tol) && NearEq(x_, other.x_, tol) &&
         NearEq(y_, other.y_, tol) && NearEq(h_, other.h_, tol);
}
template <> bool NearEq(const double& v0, const double& v1, double tol) {
  if (std::abs(v0 - v1) > tol) {
    std::cerr << "Values differ, " << v0 << " != " << v1 << std::endl;
    return false;
  }
  return true;
}

template <> bool NearEq(const Sliceable& v0, const Sliceable& v1, double tol) {
  return v0.NearEqual(v1, tol);
}

template <>
bool NearEq(const ZonoSliceInfo& v0, const ZonoSliceInfo& v1, double tol) {
  return VecsNearEq(v0.slc_vals_, v1.slc_vals_, tol);
}
template <> bool NearEq(const Interval& v0, const Interval& v1, double tol) {
  // TODO remove single &, just helpful for debugging
  return NearEq(v0.Min(), v1.Min(), tol) & NearEq(v0.Max(), v1.Max(), tol);
}

template <> bool NearEq(const Vehrs& v0, const Vehrs& v1, double tol) {
  return v0.NearEqual(v1, tol);
}

template <> bool NearEq(const FrsMega& v0, const FrsMega& v1, double tol) {
  return VecsNearEq(v0.au_, v1.au_, tol, "FrsMega::au") and
         VecsNearEq(v0.dir_, v1.dir_, tol, "FrsMega::dir") and
         VecsNearEq(v0.lan_, v1.lan_, tol, "FrsMega::lan");
}

} // namespace roahm
