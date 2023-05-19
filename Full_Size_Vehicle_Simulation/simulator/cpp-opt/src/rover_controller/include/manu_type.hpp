#ifndef ROAHM_MANU_TYPE_HPP_
#define ROAHM_MANU_TYPE_HPP_
#include <cstdint> // for uint8_t
#include <string>  // for string

/// @file manu_type.hpp Contains structures related to the maneuver type

namespace roahm {

/// Contains available maneuver types
enum class ManuType { kDirChange, kLanChange, kSpdChange, kNone };

/// Converts a string to a maneuver type, if a valid correspondence exists
/// \param manu_type_str the string to convert, should correspond with the
/// maneuver type's ToString conversion
/// \return the maneuver type
ManuType FromString(const std::string& manu_type_str);

/// Converts a maneuver type to a string, if the provided maneuver type is valid
/// \param manu_type the maneuver type to convert
/// \return a string representing the provided maneuver type
std::string ToString(const ManuType& manu_type);

// TODO this should probably be moved out of this file?
/// TODO
/// \param manu_type TODO
/// \param u0 TODO
/// \param k_param TODO
/// \return TODO
double GetGenParamAu(const ManuType& manu_type, double u0, double k_param);

// TODO this should probably be moved out of this file?
/// TODO
/// \param manu_type TODO
/// \param u0 the initial speed
/// \param k_param the Au (speed) or Ay (dir/lane change) value
/// \return TODO
double GetGenParamAy(const ManuType& manu_type, double u0, double k_param);

/// Checks whether a given maneuver type is a speed change
/// \param manu_type the maneuver type to check
/// \return true if \p manu_type is a speed change
inline bool IsSpd(ManuType manu_type) {
  return manu_type == ManuType::kSpdChange;
}

/// Checks whether a given maneuver type is a direction change
/// \param manu_type the maneuver type to check
/// \return true if \p manu_type is a direction change
inline bool IsDir(ManuType manu_type) {
  return manu_type == ManuType::kDirChange;
}

/// Checks whether a given maneuver type is a lane change
/// \param manu_type the maneuver type to check
/// \return true if \p manu_type is a lane change
inline bool IsLan(ManuType manu_type) {
  return manu_type == ManuType::kLanChange;
}

/// Checks whether a given maneuver type is a none type
/// \param manu_type the maneuver type to check
/// \return true if \p manu_type is a none type
inline bool IsNone(ManuType manu_type) { return manu_type == ManuType::kNone; }

/// Checks whether a given maneuver type is a direction or lane change
/// \param manu_type the maneuver type to check
/// \return true if \p manu_type is a direction or lane type
inline bool IsDirLan(ManuType manu_type) {
  return IsDir(manu_type) or IsLan(manu_type);
}

/// Checks whether a given maneuver type is a speed or direction change
/// \param manu_type the maneuver type to check
/// \return true if \p manu_type is a speed or direction type
inline bool IsSpdDir(ManuType manu_type) {
  return IsSpd(manu_type) or IsDir(manu_type);
}

/// Checks whether a given maneuver type is a valid maneuver type: speed, lane,
/// or direction change, not "none"
/// \param manu_type the maneuver type to check
/// \return true if \p manu_type is a valid maneuver type: a speed, lane, or
/// direction change, and not "none"
inline bool IsValid(ManuType manu_type) {
  return IsSpd(manu_type) or IsDirLan(manu_type);
}

/// Get the duration in seconds that a maneuver type is supposed to run for
/// prior to the braking portion. \param manu_type the maneuver type to check
/// \return the duration in seconds that a maneuver is supposed to run for prior
///// to the braking portion.
double GetPreBrakeTime(ManuType manu_type);

/// Returns an unsigned 8 bit integer corresponding to a maneuver type, to be
/// used in messages
/// \param manu_type the maneuver type to get an integer code for
/// \return an unsigned 8 bit integer corresponding to a maneuver type
std::uint8_t ManuToUint8(ManuType manu_type);

/// Returns a maneuver type corresponding to an integer code
/// \param manu_type the maneuver type to get an integer code for
/// \param val the value corresponding to a maneuver type
/// \return a maneuver type corresponding to the provided value
ManuType ManuFromUint8(std::uint8_t val);

} // namespace roahm

#endif // ROAHM_MANU_TYPE_HPP_
