#ifndef ROAHM_FRS_SELECT_INFO_HPP_
#define ROAHM_FRS_SELECT_INFO_HPP_

#include "manu_type.hpp"

namespace roahm {

/// The necessary information to reference a specific FRS, and whether it is
/// mirrored.
struct FrsSelectInfo {
  /// the maneuver type
  ManuType manu_type_;
  /// the u0 index number
  int idxu0_;
  /// idx0 the first index of the FRS
  int idx0_;
  /// the second index of the FRS, ignored if it is a speed change
  int idx1_;
  /// true iff the frs is mirrored
  bool mirror_;
  /// Constructor
  /// \param manu_type the maneuver type
  /// \param idxu0 the u0 index number
  /// \param idx0 the first index of the FRS
  /// \param idx1 the second index of the FRS, ignored if it is a speed change
  /// \param mirror true if the FRS is mirrored
  FrsSelectInfo(ManuType manu_type, int idxu0, int idx0, int idx1, bool mirror)
      : manu_type_{manu_type}, idxu0_{idxu0}, idx0_{idx0}, idx1_{idx1},
        mirror_{mirror} {}
};

} // namespace roahm
#endif // ROAHM_FRS_SELECT_INFO_HPP_
