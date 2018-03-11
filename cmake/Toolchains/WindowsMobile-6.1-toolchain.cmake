MESSAGE(STATUS "entering -- Windows Mobile 6 CMake Toolchain file")

# this one is important
SET(CMAKE_SYSTEM_NAME WinCE)
SET(CMAKE_SYSTEM_VERSION 5.02)
SET(CMAKE_SYSTEM_PROCESSOR ARM)

# specify the cross compiler
SET(CMAKE_C_COMPILER   "C:/Program Files/Microsoft Visual Studio 9.0/VC/ce/bin/x86_arm/cl.exe")
SET(CMAKE_CXX_COMPILER "C:/Program Files/Microsoft Visual Studio 9.0/VC/ce/bin/x86_arm/cl.exe")

# where is the target environment
SET(CMAKE_FIND_ROOT_PATH  "c:/Program Files/Windows Mobile 6 SDK/PocketPC")

# search for programs in the build host directories
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM BOTH)
# for libraries and headers in the target directories
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY BOTH)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE BOTH)

# only build in release mode for this platform
SET(CMAKE_BUILD_TYPE Release)



MESSAGE(STATUS "leaving -- Windows Mobile 6 CMake Toolchain file")