#include "risk_rtd.hpp"

#include <fmt/format.h>

#include "IpAlgTypes.hpp"
#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"
#include "fmt/core.h"
#include "frs_select_info.hpp"
#include "frs_total.hpp"
#include "ipopt_string_utils.hpp"
#include "manu_type.hpp"
#include "risk_problem_description.hpp"
#include "risk_rtd_ipopt_problem.hpp"
#include "timing_util.hpp"
#include "unit_conversion.hpp"

namespace roahm {

FrsTotal RiskRtd::GetFrsTotal() const { return frs_; }
} // namespace roahm
