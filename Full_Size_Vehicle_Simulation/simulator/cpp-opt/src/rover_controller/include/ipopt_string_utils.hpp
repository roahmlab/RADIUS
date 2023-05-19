#ifndef ROAHM_IPOPT_STRING_UTILS_HPP_
#define ROAHM_IPOPT_STRING_UTILS_HPP_

#include "IpAlgTypes.hpp"
#include "IpIpoptApplication.hpp"
#include <optional>
#include <string>

namespace roahm {

inline std::string ToString(Ipopt::SolverReturn solver_ret) {
  switch (solver_ret) {
  case Ipopt::SolverReturn::SUCCESS:
    return "SUCCESS";
  case Ipopt::SolverReturn::MAXITER_EXCEEDED:
    return "MAXITER_EXCEEDED";
  case Ipopt::SolverReturn::CPUTIME_EXCEEDED:
    return "CPUTIME_EXCEEDED";
  case Ipopt::SolverReturn::WALLTIME_EXCEEDED:
    return "WALLTIME_EXCEEDED";
  case Ipopt::SolverReturn::STOP_AT_TINY_STEP:
    return "STOP_AT_TINY_STEP";
  case Ipopt::SolverReturn::STOP_AT_ACCEPTABLE_POINT:
    return "STOP_AT_ACCEPTABLE_POINT";
  case Ipopt::SolverReturn::LOCAL_INFEASIBILITY:
    return "LOCAL_INFEASIBILITY";
  case Ipopt::SolverReturn::USER_REQUESTED_STOP:
    return "USER_REQUESTED_STOP";
  case Ipopt::SolverReturn::FEASIBLE_POINT_FOUND:
    return "FEASIBLE_POINT_FOUND";
  case Ipopt::SolverReturn::DIVERGING_ITERATES:
    return "DIVERGING_ITERATES";
  case Ipopt::SolverReturn::RESTORATION_FAILURE:
    return "RESTORATION_FAILURE";
  case Ipopt::SolverReturn::ERROR_IN_STEP_COMPUTATION:
    return "ERROR_IN_STEP_COMPUTATION";
  case Ipopt::SolverReturn::INVALID_NUMBER_DETECTED:
    return "INVALID_NUMBER_DETECTED";
  case Ipopt::SolverReturn::TOO_FEW_DEGREES_OF_FREEDOM:
    return "TOO_FEW_DEGREES_OF_FREEDOM";
  case Ipopt::SolverReturn::INVALID_OPTION:
    return "INVALID_OPTION";
  case Ipopt::SolverReturn::OUT_OF_MEMORY:
    return "OUT_OF_MEMORY";
  case Ipopt::SolverReturn::INTERNAL_ERROR:
    return "INTERNAL_ERROR";
  case Ipopt::SolverReturn::UNASSIGNED:
    return "UNASSIGNED";
  }
  return "Unknown";
}

inline std::string ToString(Ipopt::ApplicationReturnStatus status) {
  switch (status) {
  case Ipopt::ApplicationReturnStatus::Solve_Succeeded:
    return "Solve_Succeeded";
  case Ipopt::ApplicationReturnStatus::Solved_To_Acceptable_Level:
    return "Solved_To_Acceptable_Level";
  case Ipopt::ApplicationReturnStatus::Infeasible_Problem_Detected:
    return "Infeasible_Problem_Detected";
  case Ipopt::ApplicationReturnStatus::Search_Direction_Becomes_Too_Small:
    return "Search_Direction_Becomes_Too_Small";
  case Ipopt::ApplicationReturnStatus::Diverging_Iterates:
    return "Diverging_Iterates";
  case Ipopt::ApplicationReturnStatus::User_Requested_Stop:
    return "User_Requested_Stop";
  case Ipopt::ApplicationReturnStatus::Feasible_Point_Found:
    return "Feasible_Point_Found";
  case Ipopt::ApplicationReturnStatus::Maximum_Iterations_Exceeded:
    return "Maximum_Iterations_Exceeded";
  case Ipopt::ApplicationReturnStatus::Restoration_Failed:
    return "Restoration_Failed";
  case Ipopt::ApplicationReturnStatus::Error_In_Step_Computation:
    return "Error_In_Step_Computation";
  case Ipopt::ApplicationReturnStatus::Maximum_CpuTime_Exceeded:
    return "Maximum_CpuTime_Exceeded";
  case Ipopt::ApplicationReturnStatus::Maximum_WallTime_Exceeded:
    return "Maximum_WallTime_Exceeded";
  case Ipopt::ApplicationReturnStatus::Not_Enough_Degrees_Of_Freedom:
    return "Not_Enough_Degrees_Of_Freedom";
  case Ipopt::ApplicationReturnStatus::Invalid_Problem_Definition:
    return "Invalid_Problem_Definition";
  case Ipopt::ApplicationReturnStatus::Invalid_Option:
    return "Invalid_Option";
  case Ipopt::ApplicationReturnStatus::Invalid_Number_Detected:
    return "Invalid_Number_Detected";
  case Ipopt::ApplicationReturnStatus::Unrecoverable_Exception:
    return "Unrecoverable_Exception";
  case Ipopt::ApplicationReturnStatus::NonIpopt_Exception_Thrown:
    return "NonIpopt_Exception_Thrown";
  case Ipopt::ApplicationReturnStatus::Insufficient_Memory:
    return "Insufficient_Memory";
  case Ipopt::ApplicationReturnStatus::Internal_Error:
    return "Internal_Error";
  }
  return "Unknown";
};

template <typename T>
  requires requires(T x) { ToString(x); }
inline std::string ToStringOpt(const std::optional<T>& val) {
  if (val.has_value()) {
    return ToString(val.value());
  }
  return "nullopt";
}

template <typename T>
  requires requires(T x) { std::to_string(x); }
inline std::string StdToStringOpt(const std::optional<T>& val) {
  if (val.has_value()) {
    return std::to_string(val.value());
  }
  return "nullopt";
}

} // namespace roahm
#endif // ROAHM_IPOPT_STRING_UTILS_HPP_