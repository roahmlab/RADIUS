cmake_minimum_required(VERSION 3.0.2)
project(rover_control_msgs)

find_package(catkin REQUIRED COMPONENTS
  std_msgs
  geometry_msgs
  jsk_recognition_msgs
  nav_msgs 
  message_generation
)

add_message_files(
  DIRECTORY msg
  FILES RoverDebugStateStamped.msg  RoverUvrStamped.msg GenParamInfo.msg OnlineDebugMsg.msg MatlabPlotInfo.msg
)

generate_messages(
  DEPENDENCIES std_msgs jsk_recognition_msgs nav_msgs
)

catkin_package(
  CATKIN_DEPENDS std_msgs message_runtime jsk_recognition_msgs
)
