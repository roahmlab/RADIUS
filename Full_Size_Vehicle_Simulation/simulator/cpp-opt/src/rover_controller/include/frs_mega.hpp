#ifndef ROAHM_FRS_MEGA_HPP_
#define ROAHM_FRS_MEGA_HPP_

#include <iostream>
#include <vector>

#include "manu_type.hpp"
#include "vehrs.hpp"

namespace roahm {

/// TODO
struct FrsMega {
  /// The set of speed change maneuvers
  std::vector<Vehrs> au_;

  /// The set of dir change maneuvers
  std::vector<std::vector<Vehrs>> dir_;

  /// The set of lane change maneuvers
  std::vector<std::vector<Vehrs>> lan_;

  bool operator==(const FrsMega& rhs) const;

  void WriteToBinFileDirect(std::ofstream& file) const;
  static FrsMega ReadFromBinFileDirect(std::ifstream& file);

  /// Returns either the set of direction or lane changes
  /// \param manu_type the maneuver type to get, either direction or lane
  /// changes
  /// \return either the set of direction or lane changes
  const std::vector<std::vector<Vehrs>>& GetDirLanSet(ManuType manu_type) const;
};

} // namespace roahm
#endif // ROAHM_FRS_MEGA_HPP_