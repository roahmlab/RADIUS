#ifndef ROAHM_FRS_TOTAL_HPP_
#define ROAHM_FRS_TOTAL_HPP_

#include <fstream>
#include <ostream>
#include <vector>

#include <fmt/core.h>

#include "frs_mega.hpp"
#include "frs_select_info.hpp"
#include "vehrs.hpp"
namespace roahm {

/// TODO
struct FrsTotal {
  /// The initial speed ranges that each set of FRSes are compatible with
  std::vector<Interval> u0_intervals_;
  /// Sets of FRSes
  std::vector<FrsMega> megas_;
  /// The minimum initial speed
  double u0_min_;
  /// The maximum initial speed
  double u0_max_;
  /// Whether the set of FRSes was successfully loaded or not
  bool successful_;
  /// Whether only every other bin has au values
  bool alternating_au_;

  void WriteToBinFileDirect(std::ofstream& file) const;
  static FrsTotal ReadFromBinFileDirect(std::ifstream& file);

  bool NearEqual(const FrsTotal& rhs, double tol) const;
  bool operator==(const FrsTotal& other) const;

  /// Returns the minimum initial speed that can be used with the FRSes
  /// \return the minimum initial speed that can be used with the FRSes
  double MinU0() const;

  /// Returns the maximum initial speed that can be used with the FRSes
  /// \return the maximum initial speed that can be used with the FRSes
  double MaxU0() const;

  /// Returns the index of the FRSes that are compatible with the given u0
  /// \param u0 the initial speed to find a set of FRSes that are compatible
  /// with
  /// \param can_clamp_u whether an out of range u0 should be clamped to the
  /// nearest applicable range. Can be useful if noisy measurements exist
  /// \return the index of the set of FRSes that are compatible with the
  /// provided u0
  long SelectU0Idx(double& u0, bool can_clamp_u) const;

  /// Returns the index of the FRSes that are compatible with the given u0
  /// \param u0 the initial speed to find a set of FRSes that are compatible
  /// with
  /// \param can_clamp_u whether an out of range u0 should be clamped to the
  /// nearest applicable range. Can be useful if noisy measurements exist
  /// \return the index of the set of FRSes that are compatible with the
  /// provided u0
  long SelectU0IdxNoClamp(const double u0) const;

  /// Checks if only every other bin has Au values
  /// \return true iff only every other bin has Au values
  bool AlternatingAu() const;

  /// Return a reachable set corresponding with the provided information
  /// \param info the indices corresponding to a given FRS
  /// \return a reachable set corresponding with the provided information
  const Vehrs& GetVehrs(const FrsSelectInfo& info) const;

inline long GetNumFrses() const {
	long ret{0};
	for (const auto& mega : megas_) {
		ret += mega.au_.size();
		ret += mega.lan_.size();
		ret += mega.dir_.size();
	}
	return ret;
}

};


inline static std::vector<::roahm::FrsSelectInfo> FindAllValidFrses(
    const ::roahm::FrsTotal& frs_total,
    const ::roahm::RoverState predicted_state, const double max_spd = 32.0,
    const double max_pos_delta_spd = 100.0,
    const double max_neg_delta_spd = 14.0, const int every_n_spd_changes = 3,
    const bool check_mirrors = true) {
  // TODO REMOVE
  const auto u0 = predicted_state.u_;
  const auto v0 = predicted_state.v_;
  const auto r0 = predicted_state.r_;
  auto check_frs = [u0, v0, r0, max_spd, max_pos_delta_spd, max_neg_delta_spd](
                       const ::roahm::Vehrs& vehrs, const bool mirror) -> bool {
    // Used the mirrored v and r values just in case the centers
    // are ever non-zero
    const double mirror_mult = mirror ? -1 : 1;
    const auto v0_to_check = mirror_mult * v0;
    const auto r0_to_check = mirror_mult * r0;

    // Get center and generator values
    const auto u_cg = vehrs.GetUCenterGen();
    const auto v_cg = vehrs.GetVCenterGen();
    const auto r_cg = vehrs.GetRCenterGen();

    // Compute bounds
    const auto u_max = u_cg.first + std::abs(u_cg.second);
    const auto u_min = u_cg.first - std::abs(u_cg.second);
    const auto v_max = v_cg.first + std::abs(v_cg.second);
    const auto v_min = v_cg.first - std::abs(v_cg.second);
    const auto r_max = r_cg.first + std::abs(r_cg.second);
    const auto r_min = r_cg.first - std::abs(r_cg.second);

    // Check if each variable is valid
    const bool u_ok = (u0 >= u_min) and (u0 <= u_max);
    const bool v_ok = (v0_to_check >= v_min) and (v0_to_check <= v_max);
    const bool r_ok = (r0_to_check >= r_min) and (r0_to_check <= r_max);

    const bool init_conds_ok = u_ok and v_ok and r_ok;

    if (not init_conds_ok) {
      return false;
    }

    if (vehrs.GetManuType() == ManuType::kSpdChange) {
      if (vehrs.GetTrajParamMax() > max_spd) {
        return false;
      }

      if (u0 - vehrs.GetTrajParamMin() > max_neg_delta_spd) {
        return false;
      }

      // const double max_pos_dist_to_au{u0 > 19.5 ? 5 : 15.0};
      if (vehrs.GetTrajParamMax() - u0 > max_pos_delta_spd) {
        return false;
      }
    }

    return true;
  };

  std::vector<::roahm::FrsSelectInfo> frs_select_info;

  for (int u0_idx = 0; u0_idx < frs_total.megas_.size(); ++u0_idx) {
    const auto& mega = frs_total.megas_.at(u0_idx);
    const auto& au = mega.au_;
    const auto& dir = mega.dir_;
    const auto& lan = mega.lan_;
    {
      std::vector<::roahm::FrsSelectInfo> available_speed_changes;
      for (int idx0 = 0; idx0 < au.size(); ++idx0) {
        const auto& frs = au.at(idx0);
        if (check_frs(frs, false)) {
          available_speed_changes.push_back(::roahm::FrsSelectInfo{
              ::roahm::ManuType::kSpdChange, u0_idx, idx0, 0, false});
        }
      }
      if (not available_speed_changes.empty()) {
        frs_select_info.push_back(available_speed_changes.back());
        available_speed_changes.pop_back();
        for (int i = 0; i < available_speed_changes.size();
             i += every_n_spd_changes) {
          frs_select_info.push_back(available_speed_changes.at(i));
        }
      }
    }
    for (int idx0 = 0; idx0 < dir.size(); ++idx0) {
      const auto& dir_inner = dir.at(idx0);
      for (int idx1 = 0; idx1 < dir_inner.size(); ++idx1) {
        const auto& frs = dir_inner.at(idx1);
        if (check_frs(frs, false)) {
          frs_select_info.push_back(::roahm::FrsSelectInfo{
              ::roahm::ManuType::kDirChange, u0_idx, idx0, idx1, false});
        }
        if (check_mirrors and check_frs(frs, true)) {
          frs_select_info.push_back(::roahm::FrsSelectInfo{
              ::roahm::ManuType::kDirChange, u0_idx, idx0, idx1, true});
        }
      }
    }
    for (int idx0 = 0; idx0 < lan.size(); ++idx0) {
      const auto& lan_inner = lan.at(idx0);
      for (int idx1 = 0; idx1 < lan_inner.size(); ++idx1) {
        const auto& frs = lan_inner.at(idx1);
        if (check_frs(frs, false)) {
          frs_select_info.push_back(::roahm::FrsSelectInfo{
              ::roahm::ManuType::kLanChange, u0_idx, idx0, idx1, false});
        }
        if (check_mirrors and check_frs(frs, true)) {
          frs_select_info.push_back(::roahm::FrsSelectInfo{
              ::roahm::ManuType::kLanChange, u0_idx, idx0, idx1, true});
        }
      }
    }

    // Don't cover two u0 bins. If we've found any valid ones,
    // just exit before the next u0 index
    if (frs_select_info.size() > 0) {
      break;
    }
  }
  return frs_select_info;
}

} // namespace roahm
#endif // ROAHM_FRS_TOTAL_HPP_
