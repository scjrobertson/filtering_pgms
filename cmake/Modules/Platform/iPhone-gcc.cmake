SET (CMAKE_C_FLAGS_DEBUG_INIT            "${CMAKE_C_FLAGS_INIT} -DDEBUG=1 -ggdb")
SET (CMAKE_C_FLAGS_RELEASE_INIT          "${CMAKE_C_FLAGS_INIT} -DNDEBUG=1")
SET (CMAKE_C_FLAGS_RELWITHDEBINFO_INIT   "${CMAKE_C_FLAGS_INIT} -DNDEBUG=1 -ggdb")

SET (CMAKE_CXX_FLAGS_DEBUG_INIT          "${CMAKE_CXX_FLAGS_INIT} -DDEBUG=1 -ggdb")
SET (CMAKE_CXX_FLAGS_RELEASE_INIT        "${CMAKE_CXX_FLAGS_INIT} -DNDEBUG=1")
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "${CMAKE_CXX_FLAGS_INIT} -DNDEBUG=1 -ggdb")
