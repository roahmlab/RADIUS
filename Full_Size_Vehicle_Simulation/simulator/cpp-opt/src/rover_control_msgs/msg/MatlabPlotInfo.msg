Header header
bool was_successful
# Current state at time of previous trajectory start
float32 x0
float32 y0
float32 h0
float32 u0
float32 v0
float32 r0
# The predicted starting state of our trajectory
float32 pred_x
float32 pred_y
float32 pred_h 
float32 pred_u
float32 pred_v
float32 pred_r 
# Waypoint
float32 wp_x
float32 wp_y
float32 wp_h
# Chosen Values
int32 u0_idx
int32 idx0
int32 idx1
uint8 manu_type
float32 k_param
bool is_low
# Obstacles
jsk_recognition_msgs/PolygonArray obstacles 
nav_msgs/Path full_rob_path
