# iPhone Simulator CMake Toolchain definition

SET(CMAKE_SYSTEM_NAME      "iPhone")
SET(CMAKE_SYSTEM_VERSION   "${SDK_VERSION}")
SET(CMAKE_SYSTEM_PROCESSOR "i386")

SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

SET(MINVER "${DEPLOYMENT_VERSION}")
SET(SDKVER "${SDK_VERSION}")
SET(DEVTOOL "Developer") # Default
# SET(DEVTOOL "Xcode4")    # Experimental
SET(DEVROOT "/${DEVTOOL}/Platforms/iPhoneSimulator.platform/Developer")
SET(SDKROOT "${DEVROOT}/SDKs/iPhoneSimulator${SDKVER}.sdk")
# SET(CMAKE_OSX_SYSROOT "${SDKROOT}")
SET(CMAKE_OSX_SYSROOT "iphonesimulator${SDKVER}")

INCLUDE(CMakeForceCompiler)
CMAKE_FORCE_C_COMPILER("${DEVROOT}/usr/bin/gcc-4.2" iPhoneSimulator-gcc-4.2)
CMAKE_FORCE_CXX_COMPILER("${DEVROOT}/usr/bin/g++-4.2" iPhoneSimulator-g++-4.2)

SET(CMAKE_C_FLAGS_INIT   "-pipe -no-cpp-precomp --sysroot=${SDKROOT} -mmacosx-version-min=10.5 -std=c99")
SET(CMAKE_CXX_FLAGS_INIT "-pipe -no-cpp-precomp --sysroot=${SDKROOT} -mmacosx-version-min=10.5")

INCLUDE_DIRECTORIES(SYSTEM "${SDKROOT}/usr/include")
INCLUDE_DIRECTORIES(SYSTEM "${SDKROOT}/opt/iphone-simulator-${SDKVER}/include")
INCLUDE_DIRECTORIES(SYSTEM "${SDKROOT}/usr/local/iphone-simulator-${SDKVER}/include")

LINK_DIRECTORIES("${SDKROOT}/usr/lib")
LINK_DIRECTORIES("${SDKROOT}/opt/iphone-simulator-${SDKVER}/lib")
LINK_DIRECTORIES("${SDKROOT}/usr/local/iphone-simulator-${SDKVER}/lib")

SET(CMAKE_FIND_ROOT_PATH "${SDKROOT}" "/opt/iphone-simulator-${SDKVER}/" "/usr/local/iphone-simulator-${SDKVER}/")
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM BOTH)
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

SET(IPHONE 1)
SET(IPHONESIMULATOR 1)
SET(IPHONESIMULATOR_VERSION ${SDKVER})
SET(APPLE 1)
