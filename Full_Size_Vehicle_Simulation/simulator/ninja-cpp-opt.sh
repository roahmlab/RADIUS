export CMAKE_GENERATOR=Ninja && \
export C_INCLUDE_PATH=$(llvm-config-15 --includedir)
export CPLUS_INCLUDE_PATH=$(llvm-config-15 --includedir)
export LIBRARY_PATH=$(llvm-config-15 --libdir)
cd cpp-opt && \
source /opt/ros/noetic/setup.bash && \
CXXFLAGS="-isystem /usr/lib/llvm-15/lib/clang/15.0.0/include" catkin_make --use-ninja \
	-DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
	-DCMAKE_CXX_COMPILER=/usr/bin/clang++-15 \
	-DCMAKE_C_COMPILER=/usr/bin/clang-15
#	-DCMAKE_EXE_LINKER_FLAGS="-fuse-ld=mold" \
#CXXFLAGS="-isystem /usr/lib/llvm-15/lib/clang/15.0.0/include" catkin_make \
#CXXFLAGS="-isystem /usr/lib/llvm-15/lib/clang/15.0.0/include" catkin_make --use-ninja \
#	-DOpenMP_C_FLAGS=-fopenmp \
#	-DOpenMP_C_LIB_NAMES=omp \
#	-DOpenMP_omp_LIBRARY=/usr/lib/x86_64-linux-gnu/libomp5.so
#	-DOpenMP_CXX_FLAGS=-fopenmp \
#	-DOpenMP_CXX_LIB_NAMES=omp
#catkin clean -y && \
#catkin build
