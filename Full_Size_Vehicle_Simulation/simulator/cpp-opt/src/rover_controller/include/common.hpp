#ifndef ROAHM_COMMON_HPP_
#define ROAHM_COMMON_HPP_

#include <cstdint>

/// @file common.hpp Contains a common index type to use.

namespace roahm {
/// Type to be used for indices, especially those that may interact with
/// IPOPT, should generally correspond with Ipopt::Index.
using IndexT = std::size_t;
} // namespace roahm

#endif // ROAHM_COMMON_HPP_
