#ifndef ROAHM_VEHRS_HPP_
#define ROAHM_VEHRS_HPP_

#include <utility>
#include <vector>

#include "common.hpp"
#include "cuda_info.hpp"
#include "individual_mu_sigma.hpp"
#include "manu_type.hpp"
#include "mu_sigma_multi.hpp"
#include "point_xy.hpp"
#include "point_xyh.hpp"
#include "rover_state.hpp"

namespace roahm {

/// Contains information about one dimension of slicing
struct Sliceable {
  /// the dimension that the slice value is along
  int dim_;
  /// the center of the parameter's slicing range
  double center_val_;
  /// the "radius" of the slicing operation
  double slc_val_;
  /// the corresponding value in the x dimension
  double x_;
  /// the corresponding value in the y dimension
  double y_;
  /// the corresponding value in the h dimension
  double h_;
  /// Constructor
  /// \param dim The dimension to slice along
  /// \param center_val the center of the parameter's slicing range
  /// \param slc_val the "radius" of the slicing operation
  /// \param x the corresponding value in the x dimension
  /// \param y the corresponding value in the y dimension
  /// \param h the corresponding value in the h dimension
  Sliceable(int dim, double center_val, double slc_val, double x, double y,
            double h)
      : dim_(dim), center_val_(center_val), slc_val_(slc_val), x_(x), y_(y),
        h_(h) {}

  friend std::ostream& operator<<(std::ostream& out, const Sliceable& val);

  static Sliceable ReadFromBinFileDirect(std::ifstream& file);

  void WriteToBinFileDirect(std::ofstream& file) const;

  bool operator==(const Sliceable& other) const;

  bool NearEqual(const Sliceable& other, double tol) const;
};

struct SlicedInfo {
  double slc_x_;
  double slc_y_;
  double slc_h_;
  double slc_val_;
  double center_slc_val_;
  /// the sum along the x dimension from slicing
  double x_sliced_sum_;
  /// the sum along the y dimension from slicing
  double y_sliced_sum_;
  /// the sum along the h dimension from slicing
  double h_sliced_sum_;
  /// true iff the lambda values are all in \f$ [-1, 1] \f$
  bool lambda_ok_;
  /// whether the x and y values were sliced
  bool slc_xy_set_;
  SlicedInfo()
      : slc_x_{0}, slc_y_{0}, slc_h_{0}, slc_val_{0}, center_slc_val_{0},
        x_sliced_sum_{0}, y_sliced_sum_{0}, h_sliced_sum_{0}, lambda_ok_{true},
        slc_xy_set_{false} {}
};

struct ZonoSliceInfo {
  /// Output information after a slicing operation
  std::vector<Sliceable> slc_vals_;
  /// Computes the xyh values from a slicing operation
  /// \param u the longitudinal vehicle speed [m/s] to slice at
  /// \param v the lateral vehicle speed [m/s] to slice at
  /// \param r the vehicle yaw rate [rad/s] to slice at
  /// \return the computed values from the slicing operation
  SlicedInfo Slice(double u, double v, double r) const;
  bool operator==(const ZonoSliceInfo& oth) const;
};

/// A single reachable set
struct Vehrs {
  /// The xy center points of each zonotope
  std::vector<double> xy_centers_;
  /// The heading center point of each zonotope
  std::vector<double> h_centers_;
  /// The xy portion of each zonotope's generators
  std::vector<double> zono_xy_;
  /// The h portion of each zonotope's generators
  std::vector<double> zono_h_;
  /// The number of generators corresponding to each zonotope
  std::vector<IndexT> zono_sizes_;
  /// The slice values associated with each zonotope
  std::vector<ZonoSliceInfo> slc_infos_;
  /// The time intervals associated with each zonotope
  std::vector<Interval> zono_time_intervals_;
  /// Which time segment this FRS corresponds to (e.g. first or second half
  /// of a lane change)
  int t_eval_idx_;
  /// The trajectory parameter "radius" about the central value
  double k_rng_;

  std::vector<double> slc_x_;
  std::vector<double> slc_y_;
  std::vector<bool> slc_xy_set_;

  CudaInfo cuda_info_;

  [[nodiscard]] ManuType GetManuType() const noexcept(false);

  /// Associates the correct mu sigma with each zonotope's time interval
  std::vector<double> ExtractMuSigma(const MuSigmaMulti& mu_sigmas,
                                     const PointXYH& local_frame,
                                     const bool mirror) const;

  bool NearEqual(const Vehrs& oth, double tol) const;

  bool operator==(const Vehrs& oth) const;

  void WriteToBinFileDirect(std::ofstream& file) const;
  static Vehrs ReadFromBinFileDirect(std::ifstream& file);
  /// Gets the number of zonotopes in the FRS
  /// \return the number of zonotopes
  [[nodiscard]] IndexT GetNumZonos() const noexcept(true);

  Vehrs SliceAt(double u, double v, double r) const;
  Vehrs SliceAtParam(const double u, const double v, const double r,
                     const double traj_param) const;

  [[nodiscard]] double GetTrajParamMin() const noexcept(false);
  [[nodiscard]] double GetTrajParamMax() const noexcept(false);

  [[nodiscard]] std::pair<double, double> GetCenterGenTrajParam() const
      noexcept(false);

  [[nodiscard]] double GetCenterK() const noexcept(false);
  /// @brief Returns the center and generator value associated with
  /// a sliceable dimension
  /// @param target_slc_dim the sliceable dimension to look at
  /// @return the pair of center and generator values on the dimension specified
  /// by @p target_slc_dim throws if the value in that dimension is not found
  std::pair<double, double> GetCenterGen(int target_slc_dim) const
      noexcept(false);
  std::pair<double, double> GetUCenterGen() const;
  std::pair<double, double> GetVCenterGen() const;
  std::pair<double, double> GetRCenterGen() const;

  std::array<double, 3>
  GetU0V0R0SliceBeta(const ::roahm::RoverState& state) const;
  std::pair<::roahm::PointXYH, ::roahm::PointXYH>
  GetSliceableCenterAndGensNearstToTime(const double t_near) const
      noexcept(false);

  [[nodiscard]] PointXYH LerpCenter(const double t) const;

  [[nodiscard]] std::int64_t FindNearestToTimeIdx(const double t_near) const;

  [[nodiscard]] TimeMicroseconds TrajectoryDuration() const;
};
} // namespace roahm

#endif
