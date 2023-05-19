#include "dyn_obs.hpp"
#include "frs_loader.hpp"
#include "manu_type.hpp"
#include "mu_sigma_multi.hpp"
#include "risk_pdf.hpp"
#include "risk_problem_description.hpp"
#include "risk_rtd.hpp"
#include "ros/forwards.h"
#include "ros/init.h"
#include "ros/node_handle.h"
#include "ros/publisher.h"
#include "rover_state.hpp"
#include "waypoint.hpp"
#include <ros/ros.h>
#include <rover_control_msgs/GenParamInfo.h>
#include <rover_control_msgs/OnlineDebugMsg.h>
#include <stdexcept>

namespace roahm {
namespace risk_rtd_ros {

	/*
auto GenerateNewParameter(
    RiskRtd& risk_rtd,
    const RoverState& predicted_ego_state,
    const MuSigmaMulti& obs_pdf_global,
    const WaypointGlobalNoMirror& waypoint_global) {

  RiskProblemDescription risk_inputs{};
  risk_inputs.waypoint_global_no_mirror_ = waypoint_global;
  risk_inputs.rover_state_ = predicted_ego_state;
  risk_inputs.always_risky_obs_.push_back(obs_pdf_global);
  constexpr double kSmall{0.01};

  {
    const double c_x0{0.0};
    const double c_y0{-1.0 - (kSmall / 2.0)};
    const double vel{0.0};
    const double heading{0.0};
    const double len{100.0};
    const double wid{kSmall};
    const auto obs_type{DynObs::ObstacleType::kStaticBoundary};
    const DynObs static_boundary_global{
      c_x0, c_y0, heading, vel, len, wid, obs_type
    };
    risk_inputs.always_non_risky_obs_.PushObs(static_boundary_global);
  }
  return risk_rtd.RunPlanningIteration<::roahm::risk_rtd_ipopt_problem::RiskRtdIpoptProblem>(risk_inputs);
}
*/

/// Returns a ROS parameter, with a default if it is not found
/// \tparam T the output type of the parameter
/// \param name the ROS parameter name
/// \param default_val the default value if the value could not be found
/// \return the ROS parameter associated with the provided name if found,
/// otherwise returns the default
template <typename T>
T GetRosParam(const std::string name, T default_val) {
  T ret_val;
  if (not ros::param::get(name, ret_val)) {
    if constexpr (has_stream<std::ostream, T>::value) {
      ROS_ERROR_STREAM("Parameter " << name << " not found. Setting to default "
                                    << default_val);
    } else {
      ROS_ERROR_STREAM("Parameter " << name
                                    << " not found. Setting to default.");
    }
    ret_val = default_val;
  }
  return ret_val;
}


/*
class RiskRtdRos {
    const RiskRtd risk_rtd_;
    Risk_PDF risk_pdfs_;

    RiskRtdRos(const std::string frs_input_file_name) : 
      risk_rtd_{frs_input_file_name},
      risk_pdfs_{} { 
      }

    void CallbackTrajectoryStarted() {
        // The time that the next trajectory should start
        //const ros::Time predicted_next_traj_start_time_global_ = TODO;
        //// The duration of time from the start of a given scenario to the next 
        //// trajectory's starting time
        //const ros::Duration time_from_scenario_start_to_traj_start = predicted_next_traj_start_time_global_ - scenario_start_time_global_;
    }

    void CallbackInitialPlan() {
      // Callback when we want to publish a plan but the rover has not started a trajectory yet
    }

    void CallbackRoverStartedTrajectory() {
      // Callback when ego rover publishes that it has begun a new trajectory
    }


    void GenerateNewParameter() {
      RiskProblemDescription risk_inputs;
      int scenario_id_ = GetRosParam<int>("/scenario_id", -1);
      if (scenario_id_ < 0) {
        throw std::runtime_error("Scenario ID is negative!");
      }

      const std::string fname = "asdf";
      risk_pdfs_ = {};
      risk_pdfs_.GetStoredPDFs(fname, scenario_id_);


      // TODO generate world bounds here

      const RoverState obstacle_rover_state = ;
      const RoverState predicted_rover_state = PolynomialPredict();
      const double trajectory_duration = 13.0;
      const double predicted_scenario_clock_time = ;
      const auto obs_rovers_mu_sigmas = risk_pdfs_.GetMuSigmas(obstacle_rover_state, predicted_scenario_clock_time, pdf_duration);
      const auto waypoint_global = SelectWaypointGlobal(predicted_rover_state);
      const auto waypoint_local_no_mirror = waypoint_global.RelativeTo(predicted_rover_state);
      risk_inputs.always_non_risky_obs_ = {};
      risk_inputs.maybe_risky_obs_ = {};
      risk_inputs.always_risky_obs_ = {};
      risk_inputs.always_risky_obs_.emplace_back(obs_rovers_mu_sigmas);
      risk_inputs.waypoint_local_no_mirror_ = waypoint_local_no_mirror;
      risk_inputs.rover_state_ = predicted_rover_state;
      const auto out = risk_rtd_.RunPlanningIteration<::roahm::risk_rtd_ipopt_problem::RiskRtdIpoptProblem>(risk_inputs);
      PublishInputs(out);
    }
};
*/

}
}
// Params
//  /scenario_id         (int)
//  /scenario_running    (bool)
//  /scenario_start_time (double)
//  /scenario_clock_time (double) : current - start

int main(int argc, char** argv) {
	/*
    ros::init(argc, argv, "risk_rtd_ros");
    ros::NodeHandle nh;

    rover_control_msgs::GenParamInfo info_ret;
    const double u0 = 0.5;
    const double param_val = 0.3;
    const roahm::ManuType manu_type{roahm::ManuType::kSpdChange};
    const bool is_spd_change = ::roahm::IsSpd(manu_type);

    info_ret.au = is_spd_change ? param_val : u0;
    info_ret.ay = is_spd_change ? 0.0 : param_val;
    info_ret.t0_offset = 0;
    info_ret.manu_type = roahm::ManuToUint8(manu_type);
    //info_ret.header.stamp = ros::Time::now();
    ros::Publisher publisher = nh.advertise<rover_control_msgs::GenParamInfo>("/rover_7/gen_param_out", 1);

    // OK
    const int scenario_id = ::roahm::risk_rtd_ros::GetRosParam<int>("/scenario_id", -1);
    const bool scenario_running = ::roahm::risk_rtd_ros::GetRosParam<bool>("/scenario_running", false);
    const double scenario_start_time = ::roahm::risk_rtd_ros::GetRosParam<double>("/scenario_start_time", -1.0);
    const double scenario_clock_time = ::roahm::risk_rtd_ros::GetRosParam<double>("/scenario_clock_time", -1.0);
    std::cout << "Scenario Id: " << scenario_id << std::endl;
    std::cout << "Scenario Running: " << scenario_running << std::endl;
    std::cout << "Scenario Start Time: " << scenario_start_time << std::endl;
    std::cout << "Scenario Clock Time: " << scenario_clock_time << std::endl;
    while (true) {
    ros::spinOnce();
    publisher.publish(info_ret);
    }
    */
}
