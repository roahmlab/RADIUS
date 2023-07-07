cmake_minimum_required(VERSION 3.0.2)
#project(rover_controller CXX)
project(rover_controller CXX CUDA)

set(CXX_COMPILER /usr/bin/clang-15)
set(CMAKE_CXX_COMPILER /usr/bin/clang++-15)
set(CMAKE_C_COMPILER /usr/bin/clang-15)
set(CMP0079 NEW)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

#set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE include-what-you-use;--no_fwd_decls)

enable_testing()
include(CTest)
link_directories(/usr/local/lib/)

include(FetchContent)
FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG tags/9.1.0
)
FetchContent_MakeAvailable(fmt)

# 
# Google Test
#

include(FetchContent)
# FetchContent_Declare(
#   googletest
#   URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
# )
# # For Windows: Prevent overriding the parent project's compiler/linker settings
# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# FetchContent_MakeAvailable(googletest)


# Make sure that CUDA
find_package(CUDA REQUIRED)

# Make sure that MATLAB is available for MEX compilation
find_package(Matlab REQUIRED) 

# Eigen
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

# OpenMP
include(CheckLanguage)
check_language(OpenMP)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
    set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
    set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
    set(OpenMP_omp_LIBRARY /usr/local/opt/libomp/lib/libomp.dylib)
    set(OpenMP_C_INCLUDE_DIRS /usr/lib/llvm-15/include)
    set(OpenMP_CXX_INCLUDE_DIRS /usr/lib/llvm-15/include)
    set(OPENMP_FOUND true)
  endif()
endif()
#find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  message("OPENMP_FOUND")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wall")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Catkin for ROS integration
find_package(catkin REQUIRED
        roscpp message_generation message_runtime std_msgs geometry_msgs joy tf2_msgs rover_control_msgs)



set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_DOXYGEN_DOCS ON)


if (BUILD_DOXYGEN_DOCS)
  find_package(Doxygen REQUIRED)
endif()

set(HOME $ENV{HOME})
set(CMAKE_CXX_STANDARD 20)
set(CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
set(CONTROLLER_NAME ${CONTROLLER_NODE_NAME})
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake")


set(IPOPT_INCLUDE_DIRS "/usr/local/include/coin-or/")


#add_library(ipopt STATIC IMPORTED)
#set_target_properties(i)
#find_library(ipopt NAMES libipopt.a PATHS /usr/local/lib)
#if (IPOPT_FOUND)
#  message("FOUND IPOPT")
#else()
#  message("DIDNT FIND IPOPT")
#endif()


set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_EXE_LINKER_FLAGS_OLD "${CMAKE_EXE_LINKER_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "-static")
find_library(ipopt STATIC IMPORTED)
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
set(BUILD_SHARED_LIBS ON)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS_OLD}")

include_directories(${IPOPT_INCLUDE_DIRS})

catkin_package(
  INCLUDE_DIRS ${IPOPT_INCLUDE_DIRS}
  LIBRARIES rover_controller
  CATKIN_DEPENDS std_msgs geometry_msgs tf2_msgs)

# For files
file(GLOB_RECURSE ALL_SOURCE_FILES "src/*.cc" "include/*.hpp")
file(GLOB_RECURSE INCLUDE_FILES "include/*.hpp")
file(GLOB_RECURSE SRC_FILES_NO_EXEC "src/*.cc")

add_library(CudaLib cuda_src/prob_integ.cu cuda_src/pre_slice.cu)
target_include_directories(CudaLib PUBLIC cuda_include include ${CUDA_INCLUDE_DIRS})
target_link_libraries(CudaLib ${CMAKE_CUDA_RUNTIME_LIBRARY})
target_compile_options(CudaLib PUBLIC -O3)

add_library(MainLib ${SRC_FILES_NO_EXEC})
target_include_directories(MainLib PUBLIC src include ${catkin_INCLUDE_DIRS} /usr/include/python3.8 ${CUDA_INCLUDE_DIRS})
target_link_libraries(MainLib CudaLib ${catkin_LIBRARIES} ipopt python3.8 pthread Eigen3::Eigen fmt::fmt)

function(add_rover_executable EXEC_TARGET_NAME EXEC_FNAME IS_MEX)
  if(IS_MEX)
    matlab_add_mex(NAME ${EXEC_TARGET_NAME} SRC ${EXEC_FNAME})
  else()
    add_executable(${EXEC_TARGET_NAME} ${EXEC_FNAME})
  endif()

  add_dependencies(${EXEC_TARGET_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

  target_compile_options(${EXEC_TARGET_NAME} PUBLIC -O3 -fopenmp)
  target_include_directories(${EXEC_TARGET_NAME} PUBLIC ${catkin_INCLUDE_DIRS} ${IPOPT_INCLUDE_DIRS} /usr/include/python3.8 )
  target_link_libraries(${EXEC_TARGET_NAME} ${catkin_LIBRARIES} ipopt python3.8 pthread MainLib Eigen3::Eigen)
endfunction()

function(add_rover_test EXEC_TARGET_NAME EXEC_FNAME)
  add_executable(${EXEC_TARGET_NAME} ${EXEC_FNAME} ${INCLUDE_FILES})
  add_dependencies(${EXEC_TARGET_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
  target_include_directories(${EXEC_TARGET_NAME} PUBLIC ${catkin_INCLUDE_DIRS} ${IPOPT_INCLUDE_DIRS} /usr/include/python3.8 )
  target_link_libraries(${EXEC_TARGET_NAME} MainLib)
  #target_link_libraries(${EXEC_TARGET_NAME} ${catkin_LIBRARIES} ipopt python3.8 pthread MainLib Eigen3::Eigen)
  target_compile_options(${EXEC_TARGET_NAME} PUBLIC -O3 -fopenmp -Wall )
  target_link_libraries(${EXEC_TARGET_NAME} gtest gtest_main)
  include(GoogleTest)
  gtest_discover_tests(${EXEC_TARGET_NAME})
endfunction()

#message("CUDA INCLUDE: ${CUDA_INCLUDE_DIRS}")
#message("CUDA RUNTIME LIBRARIES: ${CUDA_RUNTIME_LIBRARIES}")
#message("CUDA LIBRARIES: ${CUDA_LIBRARIES}}")

add_rover_test(GENCON_TEST test/gencon_test.cc)
add_rover_test(RISK_RTD_IPOPT_TEST test/risk_rtd_ipopt_test.cc)
add_rover_test(DYN_OBS_TEST test/dyn_obs_test.cc)
add_rover_test(FRS_IO_TEST test/frs_io_test.cc)
add_rover_test(PROCESS_UTILS_TEST test/process_utils_test.cc)
add_rover_test(CUDA_PROB_TEST test/cuda_prob_test.cc)
add_rover_test(FL_ZONO_CONSTRAINT_TEST test/fl_zono_constraint_test.cc)
add_rover_test(INDIVIDUAL_MU_SIGMA_TEST test/individual_mu_sigma_test.cc)
add_rover_test(WAYPOINT_COST_TEST test/waypoint_cost_test.cc)
add_rover_executable(RISK_RTD_MEX mex/risk_rtd_mex.cc TRUE)
add_rover_executable(MATLAB_PROCESS_FRS mex/matlab_process_frses.cc TRUE)
add_rover_executable(RISK_RTD_ROS ros_exec/risk_rtd_ros.cc FALSE)

#if (BUILD_DOXYGEN_DOCS)
#  add_subdirectory(docs)
#endif()
