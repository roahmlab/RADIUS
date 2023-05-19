#include "vehrs.hpp"
#include <algorithm>
#include <fmt/format.h>

#include <stdexcept>

#include "frs_io.hpp"
#include "manu_type.hpp"
#include "mu_sigma_multi.hpp"
namespace roahm {
namespace {
ManuType
CudaInfoManuTypeToManuType(const int cuda_info_manu_type) noexcept(false) {
  // TODO just move this into the loading phase
  fmt::print("CudaInfoManuTypeToManuType: {}\n", cuda_info_manu_type);
  if (cuda_info_manu_type == 1) {
    return ManuType::kSpdChange;
  } else if (cuda_info_manu_type == 2) {
    return ManuType::kDirChange;
  } else if (cuda_info_manu_type == 3) {
    return ManuType::kLanChange;
  }
  throw std::runtime_error(
      "Unable to determine manuever type of cuda info manu type value " +
      std::to_string(cuda_info_manu_type));
}
} // namespace
std::ostream& operator<<(std::ostream& out, const Sliceable& val) {
  return out << "\n"
             << "x: " << val.x_ << "\n"
             << "y: " << val.y_ << "\n"
             << "h: " << val.h_ << "\n"
             << "center_val_: " << val.center_val_ << "\n"
             << "slc_val_: " << val.slc_val_ << "\n"
             << "dim: " << val.dim_ << "\n";
}
Sliceable Sliceable::ReadFromBinFileDirect(std::ifstream& file) {
  // TODO move to impl, use proper fcns
  int dim;
  double center_val, slc_val, x, y, h;
  file.read(reinterpret_cast<char*>(&dim), sizeof(int));
  file.read(reinterpret_cast<char*>(&center_val), sizeof(double));
  file.read(reinterpret_cast<char*>(&slc_val), sizeof(double));
  file.read(reinterpret_cast<char*>(&x), sizeof(double));
  file.read(reinterpret_cast<char*>(&y), sizeof(double));
  file.read(reinterpret_cast<char*>(&h), sizeof(double));
  return Sliceable(dim, center_val, slc_val, x, y, h);
}
void Sliceable::WriteToBinFileDirect(std::ofstream& file) const {
  // TODO move to impl, use proper fcns
  file.write(reinterpret_cast<const char*>(&dim_), sizeof(int));
  file.write(reinterpret_cast<const char*>(&center_val_), sizeof(double));
  file.write(reinterpret_cast<const char*>(&slc_val_), sizeof(double));
  file.write(reinterpret_cast<const char*>(&x_), sizeof(double));
  file.write(reinterpret_cast<const char*>(&y_), sizeof(double));
  file.write(reinterpret_cast<const char*>(&h_), sizeof(double));
}
bool Sliceable::operator==(const Sliceable& other) const {
  return dim_ == other.dim_ && center_val_ == other.center_val_ &&
         slc_val_ == other.slc_val_ && x_ == other.x_ && y_ == other.y_ &&
         h_ == other.h_;
}

// TODO we should just make this a member variable and handle at loading time
[[nodiscard]] ManuType Vehrs::GetManuType() const noexcept(false) {
  return CudaInfoManuTypeToManuType(cuda_info_.manu_type_);
}

[[nodiscard]] std::pair<double, double> Vehrs::GetCenterGenTrajParam() const
    noexcept(false) {
  return {cuda_info_.cg_p_.at(0), cuda_info_.cg_p_.at(1)};
}
[[nodiscard]] double Vehrs::GetTrajParamMin() const noexcept(false) {
  return cuda_info_.cg_p_.at(0) - std::abs(cuda_info_.cg_p_.at(1));
}
[[nodiscard]] double Vehrs::GetTrajParamMax() const noexcept(false) {
  return cuda_info_.cg_p_.at(0) + std::abs(cuda_info_.cg_p_.at(1));
}
[[nodiscard]] IndexT Vehrs::GetNumZonos() const noexcept(true) {
  return zono_sizes_.size();
}

std::vector<double> Vehrs::ExtractMuSigma(const MuSigmaMulti& mu_sigmas,
                                          const PointXYH& local_frame,
                                          const bool mirror = false) const {
  const int expected_out_size =
      MuSigmaMulti::kEltsPerMuSigma * zono_time_intervals_.size();

  const int min_num_mu_sigmas = mu_sigmas.MinNumMuSigmas();
  if (mu_sigmas.MinNumMuSigmas() * MuSigmaMulti::kEltsPerMuSigma <
      expected_out_size) {
    throw std::runtime_error(
        "Mu sigma does not contain enough elements to cover the FRS: " +
        std::to_string(min_num_mu_sigmas) + " < " +
        std::to_string(expected_out_size));
  }

  std::vector<double> ret;
  ret.reserve(expected_out_size);
  for (int i = 0; i < zono_time_intervals_.size(); i++) {
    const TimeMicroseconds curr_width_us{
        TimeMicroseconds::FromSeconds(zono_time_intervals_.at(i).Width())};
    const auto& mu_sigma_to_pull_from = mu_sigmas.GetMuSigmaVec(curr_width_us);
    const auto curr_interval_min_us =
        TimeMicroseconds::FromSeconds(zono_time_intervals_.at(i).Min());
    const auto us_between_mu_sigmas = mu_sigmas.MicrosecondsBetweenMuSigmas();
    const auto vec_idx = curr_interval_min_us / us_between_mu_sigmas;

    const int start_idx = vec_idx * MuSigmaMulti::kEltsPerMuSigma;
    const double mu_x = mu_sigma_to_pull_from.at(start_idx + 0);
    const double mu_y = mu_sigma_to_pull_from.at(start_idx + 1);
    const double sigma_1 = mu_sigma_to_pull_from.at(start_idx + 2);
    const double sigma_2 = mu_sigma_to_pull_from.at(start_idx + 3);
    const double sigma_4 = mu_sigma_to_pull_from.at(start_idx + 4);
    const ::roahm::IndividualMuSigma mu_sigma_global{mu_x, mu_y, sigma_1,
                                                     sigma_2, sigma_4};
    const auto mu_sigma_local = mu_sigma_global.RelativeTo(local_frame, mirror);

    ret.push_back(mu_sigma_local.mu_x_);
    ret.push_back(mu_sigma_local.mu_y_);
    ret.push_back(mu_sigma_local.sigma_1_);
    ret.push_back(mu_sigma_local.sigma_2_);
    ret.push_back(mu_sigma_local.sigma_4_);
  }

  if (ret.size() != expected_out_size) {
    throw std::runtime_error(
        "Expected mu sigma size: " + std::to_string(expected_out_size) +
        " but got: " + std::to_string(ret.size()));
  }
  return ret;
}

bool Vehrs::NearEqual(const Vehrs& oth, double tol) const {
  const bool t_eval_idx_eq = t_eval_idx_ == oth.t_eval_idx_;
  if (not t_eval_idx_eq) {
    std::cout << "t_eval_idx_ not equal"
              << "[" << t_eval_idx_ << " != " << oth.t_eval_idx_ << "]"
              << std::endl;
  }
  constexpr double kTrajParamTol = 1e-4;
  const bool k_rng_near_eq = NearEq(k_rng_, oth.k_rng_, kTrajParamTol);
  if (not k_rng_near_eq) {
    std::cout << "k_rng_ not near equal"
              << "[" << k_rng_ << " != " << oth.k_rng_ << "]" << std::endl;
  }
  return VecsNearEq(xy_centers_, oth.xy_centers_, tol,
                    "Vehrs::xy_centers_ [tol]") and
         VecsNearEq(h_centers_, oth.h_centers_, tol,
                    "Vehrs::h_centers_ [tol]") and
         VecsNearEq(zono_xy_, oth.zono_xy_, tol, "Vehrs::zono_xy_ [tol]") and
         VecsNearEq(zono_h_, oth.zono_h_, tol, "Vehrs::zono_h [tol]") and
         VecsEq(zono_sizes_, oth.zono_sizes_, "Vehrs::zono_sizes_ [==]") and
         VecsNearEq(slc_infos_, oth.slc_infos_, tol,
                    "Vehrs::slc_infos_ [tol]") and
         VecsNearEq(zono_time_intervals_, oth.zono_time_intervals_, tol,
                    "Vehrs::zono_time_intervals_ [tol]") and
         t_eval_idx_eq and k_rng_near_eq and
         VecsNearEq(slc_x_, oth.slc_x_, tol, "Vehrs::slc_x_ [tol]") and
         VecsNearEq(slc_y_, oth.slc_y_, tol, "Vehrs::slc_y_ [tol]") and
         VecsEq(slc_xy_set_, oth.slc_xy_set_) and
         (cuda_info_ == oth.cuda_info_);
}
bool Vehrs::operator==(const Vehrs& oth) const {
  return VecsEq(xy_centers_, oth.xy_centers_, "Vehrs::xy_centers") and
         VecsEq(h_centers_, oth.h_centers_, "Vehrs::h_centers") and
         VecsEq(zono_xy_, oth.zono_xy_, "Vehrs::zono_xy") and
         VecsEq(zono_h_, oth.zono_h_, "Vehrs::zono_h") and
         VecsEq(zono_sizes_, oth.zono_sizes_, "Vehrs::zono_sizes") and
         VecsEq(slc_infos_, oth.slc_infos_, "Vehrs::slc_infos") and
         VecsEq(zono_time_intervals_, oth.zono_time_intervals_,
                "Vehrs::zono_time_intervals") and
         (t_eval_idx_ == oth.t_eval_idx_) and (k_rng_ == oth.k_rng_) and
         VecsEq(slc_x_, oth.slc_x_, "Vehrs::slc_x") and
         VecsEq(slc_y_, oth.slc_y_, "Vehrs::slc_y") and
         VecsEq(slc_xy_set_, oth.slc_xy_set_, "Vehrs::slc_xy_set") and
         (cuda_info_ == oth.cuda_info_);
}

std::pair<double, double> Vehrs::GetCenterGen(int target_slc_dim) const
    noexcept(false) {
  const auto& slc_info = slc_infos_.at(0);
  const auto& slc_vals = slc_info.slc_vals_;
  for (const auto& sliceable : slc_vals) {
    const auto sl_dim = sliceable.dim_;
    const auto center_val = sliceable.center_val_;
    const auto gen_val = sliceable.slc_val_;
    if (sl_dim == target_slc_dim) {
      return {center_val, gen_val};
    }
  }
  throw std::runtime_error("No center/generator value found on dimension " +
                           std::to_string(target_slc_dim));
  return {0.0, 0.0};
}
[[nodiscard]] double Vehrs::GetCenterK() const noexcept(false) {
  const auto& slc_info = slc_infos_.at(0);
  const auto& slc_vals = slc_info.slc_vals_;
  for (const auto& sliceable : slc_vals) {
    const auto sl_dim = sliceable.dim_;
    const auto center_val = sliceable.center_val_;
    if (sl_dim == 11 || sl_dim == 12) {
      return center_val;
    }
  }
  throw std::runtime_error("No center trajectory parameter value found");
  return 0.0;
}

std::pair<double, double> Vehrs::GetUCenterGen() const {
  return GetCenterGen(7);
}
std::pair<double, double> Vehrs::GetVCenterGen() const {
  return GetCenterGen(8);
}
std::pair<double, double> Vehrs::GetRCenterGen() const {
  return GetCenterGen(9);
}

std::array<double, 3>
Vehrs::GetU0V0R0SliceBeta(const ::roahm::RoverState& state) const {
  const auto u_cen_gen = GetUCenterGen();
  const auto v_cen_gen = GetVCenterGen();
  const auto r_cen_gen = GetRCenterGen();
  return {(state.u_ - u_cen_gen.first) / u_cen_gen.second,
          (state.v_ - v_cen_gen.first) / v_cen_gen.second,
          (state.r_ - r_cen_gen.first) / r_cen_gen.second};
}

std::int64_t Vehrs::FindNearestToTimeIdx(const double t_near) const {
  double min_t_dist = std::numeric_limits<double>::max();
  std::int64_t min_t_idx = -1;
  for (int i = 0; i < zono_time_intervals_.size(); ++i) {
    const auto curr_dist = zono_time_intervals_.at(i).DistanceTo(t_near);
    if (curr_dist < min_t_dist) {
      min_t_dist = curr_dist;
      min_t_idx = i;
    }
  }
  return min_t_idx;
}
PointXYH Vehrs::LerpCenter(const double t) const {
  if (zono_time_intervals_.empty()) {
    return {};
  }
  if (zono_time_intervals_.size() == 1) {
    return {xy_centers_.at(0), xy_centers_.at(1), h_centers_.at(0)};
  }

  const auto min_t_idx{FindNearestToTimeIdx(t)};

  if (min_t_idx < 0 or min_t_idx >= zono_time_intervals_.size()) {
    throw std::runtime_error("Invalid min t idx");
  }

  const std::int64_t num_zonos{
      static_cast<std::int64_t>(zono_time_intervals_.size())};
  const auto t_center_nearest{zono_time_intervals_.at(min_t_idx).Midpoint()};
  const std::int64_t t_next_idx{[&]() -> std::int64_t {
    if (t - t_center_nearest > 0) {
      // t is to the right the nearest t center
      if (min_t_idx + 1 < num_zonos) {
        return min_t_idx + 1;
      }
      return min_t_idx - 1;
    } else {
      // t is to the right the nearest t center
      if (min_t_idx - 1 >= 0) {
        // We can go to the left
        return min_t_idx - 1;
      }
      return min_t_idx + 1;
    }
  }()};

  const double t_center_next{zono_time_intervals_.at(t_next_idx).Midpoint()};
  // fmt::print("t = {}, t_center_nearest = {}, t_center_next = {}, num = {},
  // denom = {}\n", t, t_center_nearest, t_center_next, t - t_center_nearest);
  const double alpha{(t - t_center_nearest) /
                     (t_center_next - t_center_nearest)};

  const double x0{xy_centers_.at(min_t_idx * 2 + 0)};
  const double y0{xy_centers_.at(min_t_idx * 2 + 1)};
  const double h0{h_centers_.at(min_t_idx)};

  const double x1{xy_centers_.at(t_next_idx * 2 + 0)};
  const double y1{xy_centers_.at(t_next_idx * 2 + 1)};
  const double h1{h_centers_.at(t_next_idx)};

  const double x_c{(alpha * x1) + ((1.0 - alpha) * x0)};
  const double y_c{(alpha * y1) + ((1.0 - alpha) * y0)};
  const double h_c{(alpha * h1) + ((1.0 - alpha) * h0)};
  // fmt::print("Lerping betwen ({}, {}, t={}) and ({}, {}, t={}) -> ({}, {},
  // alpha={})\n", x0, y0, t_center_nearest, x1, y1, t_center_next, x_c, y_c,
  // alpha);

  return {x_c, y_c, h_c};
}

std::pair<::roahm::PointXYH, ::roahm::PointXYH>
Vehrs::GetSliceableCenterAndGensNearstToTime(const double t_near) const
    noexcept(false) {
  const auto min_t_idx{FindNearestToTimeIdx(t_near)};
  const auto& zono_slc_info = slc_infos_.at(min_t_idx);
  for (const auto& slc_val : zono_slc_info.slc_vals_) {
    // TODO make this more robust
    if (slc_val.dim_ == 11 or slc_val.dim_ == 12) {
      return {PointXYH{xy_centers_.at(min_t_idx * 2),
                       xy_centers_.at(min_t_idx * 2 + 1),
                       h_centers_.at(min_t_idx)},
              PointXYH{slc_val.x_, slc_val.y_, slc_val.h_}};
    }
  }
  throw std::runtime_error("No sliceable value found");
  return {PointXYH{0.0, 0.0, 0.0}, PointXYH{0.0, 0.0, 0.0}};
}

Vehrs Vehrs::SliceAt(double u, double v, double r) const {
  Vehrs ret{};
  ret.cuda_info_ = cuda_info_;
  ret.zono_time_intervals_ = zono_time_intervals_;
  ret.xy_centers_.resize(xy_centers_.size());
  ret.h_centers_.resize(h_centers_.size());
  ret.slc_x_.resize(h_centers_.size());
  ret.slc_y_.resize(h_centers_.size());
  ret.slc_xy_set_.resize(h_centers_.size());
  for (std::size_t i = 0; i < zono_sizes_.size(); ++i) {
    const auto sliced_info = slc_infos_.at(i).Slice(u, v, r);
    const double unsliced_center_x = xy_centers_.at(i * 2);
    const double unsliced_center_y = xy_centers_.at(i * 2 + 1);
    const double zono_center_x = unsliced_center_x + sliced_info.x_sliced_sum_;
    const double zono_center_y = unsliced_center_y + sliced_info.y_sliced_sum_;
    const double zono_center_h = sliced_info.h_sliced_sum_ + h_centers_.at(i);
    ret.xy_centers_.at(i * 2) = zono_center_x;
    ret.xy_centers_.at(i * 2 + 1) = zono_center_y;
    ret.h_centers_.at(i) = zono_center_h;

    for (const auto& sl : slc_infos_.at(i).slc_vals_) {
      if (sl.dim_ == 11 || sl.dim_ == 12) {
        ret.slc_x_.at(i) = sl.x_;
        ret.slc_y_.at(i) = sl.y_;
        ret.slc_xy_set_.at(i) = true;
      }
    }
    if (not ret.slc_xy_set_.at(i)) {
      throw std::runtime_error("Not sliced");
    }
  }
  ret.zono_sizes_ = zono_sizes_;
  ret.zono_h_ = zono_h_;
  ret.zono_xy_ = zono_xy_;
  ret.slc_infos_ = slc_infos_;
  return ret;
}

Vehrs Vehrs::SliceAtParam(const double u, const double v, const double r,
                          const double traj_param) const {
  Vehrs ret{};
  ret.cuda_info_ = cuda_info_;
  ret.zono_time_intervals_ = zono_time_intervals_;
  ret.xy_centers_.resize(xy_centers_.size());
  ret.h_centers_.resize(h_centers_.size());
  ret.slc_x_.resize(h_centers_.size());
  ret.slc_y_.resize(h_centers_.size());
  ret.slc_xy_set_.resize(h_centers_.size());
  for (std::size_t i = 0; i < zono_sizes_.size(); ++i) {
    const auto sliced_info = slc_infos_.at(i).Slice(u, v, r);
    double param_x_offset{0.0};
    double param_y_offset{0.0};
    for (const auto& sl : slc_infos_.at(i).slc_vals_) {
      if (sl.dim_ == 11 || sl.dim_ == 12) {
        const double lambda = (traj_param - sl.center_val_) / sl.slc_val_;
        if (std::abs(lambda) - 1.0 >= 1.0e-3) {
          std::cout << "Traj Param: " << traj_param << std::endl;
          std::cout << "Center: " << sl.center_val_ << std::endl;
          std::cout << "SLC VAL: " << sl.slc_val_ << std::endl;
          throw std::runtime_error("Invalid slice param");
        }
        param_x_offset = lambda * sl.x_;
        param_y_offset = lambda * sl.y_;

        ret.slc_x_.at(i) = sl.x_;
        ret.slc_y_.at(i) = sl.y_;
        ret.slc_xy_set_.at(i) = true;
      }
    }
    if (not(ret.slc_xy_set_.at(i))) {
      throw std::runtime_error("Not sliced.");
    }
    const double unsliced_center_x = xy_centers_.at(i * 2);
    const double unsliced_center_y = xy_centers_.at(i * 2 + 1);
    const double x_offset = sliced_info.x_sliced_sum_ + param_x_offset;
    const double y_offset = sliced_info.y_sliced_sum_ + param_y_offset;
    const double zono_center_x = unsliced_center_x + x_offset;
    const double zono_center_y = unsliced_center_y + y_offset;
    const double zono_center_h = sliced_info.h_sliced_sum_ + h_centers_.at(i);
    ret.xy_centers_.at(i * 2) = zono_center_x;
    ret.xy_centers_.at(i * 2 + 1) = zono_center_y;
    ret.h_centers_.at(i) = zono_center_h;
  }
  ret.zono_sizes_ = zono_sizes_;
  ret.zono_h_ = zono_h_;
  ret.zono_xy_ = zono_xy_;
  ret.slc_infos_ = slc_infos_;
  return ret;
}

// Vehrs
void Vehrs::WriteToBinFileDirect(std::ofstream& file) const {
  WriteToBinFile(file, xy_centers_);
  WriteToBinFile(file, h_centers_);
  WriteToBinFile(file, zono_xy_);
  WriteToBinFile(file, zono_h_);
  WriteToBinFile(file, zono_sizes_);
  WriteToBinFile(file, slc_infos_);
  WriteToBinFile(file, zono_time_intervals_);
  WriteToBinFile(file, t_eval_idx_);
  WriteToBinFile(file, k_rng_);
  WriteToBinFile(file, slc_x_);
  WriteToBinFile(file, slc_y_);
  WriteToBinFile(file, slc_xy_set_);
  WriteToBinFile(file, cuda_info_);
}

Vehrs Vehrs::ReadFromBinFileDirect(std::ifstream& file) {
  Vehrs ret;
  ReadToVarFromBinFile(ret.xy_centers_, file);
  ReadToVarFromBinFile(ret.h_centers_, file);
  ReadToVarFromBinFile(ret.zono_xy_, file);
  ReadToVarFromBinFile(ret.zono_h_, file);
  ReadToVarFromBinFile(ret.zono_sizes_, file);
  ReadToVarFromBinFile(ret.slc_infos_, file);
  ReadToVarFromBinFile(ret.zono_time_intervals_, file);
  ReadToVarFromBinFile(ret.t_eval_idx_, file);
  ReadToVarFromBinFile(ret.k_rng_, file);
  ReadToVarFromBinFile(ret.slc_x_, file);
  ReadToVarFromBinFile(ret.slc_y_, file);
  ReadToVarFromBinFile(ret.slc_xy_set_, file);
  ReadToVarFromBinFile(ret.cuda_info_, file);
  return ret;
}

[[nodiscard]] TimeMicroseconds Vehrs::TrajectoryDuration() const {
  return ::roahm::TimeMicroseconds::FromSeconds(
      zono_time_intervals_.back().Max());
}

} // namespace roahm
