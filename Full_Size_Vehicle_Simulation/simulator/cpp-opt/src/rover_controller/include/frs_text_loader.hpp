#ifndef ROAHM_FRS_TEXT_LOADER_HPP_
#define ROAHM_FRS_TEXT_LOADER_HPP_

#include "frs_total.hpp"

namespace roahm {

/// Loads FRSes from a file
/// \param f_name the file name to load
/// \return the loaded FRS information, which contains a flag for whether the
/// operation was successful or not
[[deprecated]] FrsTotal LoadFrs(const std::string& f_name);

} // namespace roahm
#endif // ROAHM_FRS_TEXT_LOADER_HPP_