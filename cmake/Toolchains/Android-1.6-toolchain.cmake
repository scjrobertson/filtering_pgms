# Android NDK 1.6 CMake Toolchain definition

SET(CMAKE_SYSTEM_NAME      "Android")
SET(CMAKE_SYSTEM_VERSION   "1.6")
SET(CMAKE_SYSTEM_PROCESSOR "arm")

SET(NDK_ROOT         "$ENV{ANDROID_NDK}")
SET(NDKVER           "1.6")
SET(TARGET_PLATFORM  "android-4")
SET(TARGET_ARCH_ABI  "arm")
SET(TARGET_ABI       "${TARGET_PLATFORM}-${TARGET_ARCH_ABI}")

# The SYSROOT point to a directory
# that contains all public header files for a given platform, plus
# some libraries and object files used for linking the generated
# target files properly.
SET(SYSROOT          "${NDK_ROOT}/build/platforms/${TARGET_PLATFORM}/arch-${TARGET_ARCH_ABI}")

# specific toolchain
IF(APPLE)
  SET(HOST_PREBUILT    "${NDK_ROOT}/build/prebuilt/darwin-x86")
ELSEIF(LINUX)
  SET(HOST_PREBUILT    "${NDK_ROOT}/build/prebuilt/linux-x86")
ENDIF(APPLE)
SET(TOOLCHAIN_NAME   "arm-eabi-4.2.1")
SET(TOOLCHAIN_PREFIX "${HOST_PREBUILT}/${TOOLCHAIN_NAME}/bin/arm-eabi-")

# specify the cross compiler
SET(CMAKE_C_COMPILER    "${TOOLCHAIN_PREFIX}gcc")
SET(CMAKE_CXX_COMPILER  "${TOOLCHAIN_PREFIX}g++")

# include directories
INCLUDE_DIRECTORIES(SYSTEM "${SYSROOT}/usr/include")

# link directories
LINK_DIRECTORIES("${SYSROOT}/usr/lib")

# where is the target environment
SET(CMAKE_FIND_ROOT_PATH "${HOST_PREBUILT}/${TOOLCHAIN_NAME}")

# search for programs in the build host directories
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM BOTH)
# for libraries and headers in the target directories
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY BOTH)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE BOTH)

# only build in release mode for this platform
SET(CMAKE_BUILD_TYPE Release)
