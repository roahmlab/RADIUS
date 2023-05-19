cmake_minimum_required(VERSION 3.16)
find_package(Doxygen REQUIRED)

# Find all the public headers
get_target_property(ROVER_RTD_PUBLIC_HEADER_DIR RISK_RTD_MEX INTERFACE_INCLUDE_DIRECTORIES)
file(GLOB_RECURSE ROVER_RTD_PUBLIC_HEADERS ${ROVER_RTD_PUBLIC_HEADER_DIR}/*.hpp)

set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/include)
message("DOXYGEN_INPUT_DIR = ${DOXYGEN_INPUT_DIR}")
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
message("DOXYFILE_IN = ${DOXYFILE_IN}")
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

#Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)
message("DOXYFILE_OUT = ${DOXYFILE_OUT}")
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us
add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
        DEPENDS ${ROVER_RTD_PUBLIC_HEADERS}
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
        MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
        COMMENT "Generating docs")

# Nice named target so we can run the job easily
add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})