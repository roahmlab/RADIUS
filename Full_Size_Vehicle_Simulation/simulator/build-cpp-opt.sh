cd cpp-opt/src && \
source /opt/ros/noetic/setup.bash && \
catkin build -DCMAKE_C_COMPILER="/usr/bin/clang-15" -DCMAKE_CXX_COMPILER="/usr/bin/clang++-15"
