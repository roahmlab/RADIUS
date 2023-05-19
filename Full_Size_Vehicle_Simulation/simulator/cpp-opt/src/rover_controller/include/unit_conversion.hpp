#ifndef ROAHM_UNIT_CONVERSION_HPP_
#define ROAHM_UNIT_CONVERSION_HPP_
namespace roahm {
constexpr inline double
MilesPerHourToMetersPerSecond(const double miles_per_hour) {
  return miles_per_hour / 2.2369362921;
}
constexpr inline double
MetersPerSecondToMilesPerHour(const double meters_per_second) {
  return meters_per_second * 2.2369362921;
}
} // namespace roahm
#endif // ROAHM_UNIT_CONVERSION_HPP_
