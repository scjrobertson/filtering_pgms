#
# filtering_pgms/src CMakeLists.txt
#

#Enable warnings
add_definitions("-Wall")

#Include the headers and src folders
file(GLOB HEADERS "${CMAKE_SOURCE_DIR}/include/*.hpp")
file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cc")

message(${HEADERS})

#Create a library for unit testing
add_library(filters_base ${HEADERS} ${SOURCES})

#Compile the code
add_executable(filters main.cc)
target_link_libraries(filters filters_base emdw ${GLINEAR_LIBRARIES} ${PATRECII_LIBRARIES} ${BOOST_LIBRARIES})
