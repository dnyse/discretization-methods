cmake_minimum_required(VERSION 3.16)
project(solver CXX)
include(FetchContent)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_compile_definitions(_Alignof=__alignof)

# Set different flags for different build types
# if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
#     message(STATUS "Configuring Debug build")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -Wextra")
# elseif(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
#     message(STATUS "Configuring Release build")
# endif()

# Check for Boost 1.87.0 as a required dependency
message(STATUS "Checking for Boost 1.87.0 installation...")

# First try to find Boost using CMake's config mode (preferred method)
set(Boost_NO_SYSTEM_PATHS OFF)
set(Boost_NO_BOOST_CMAKE OFF)
find_package(Boost 1.87.0 CONFIG)

if(NOT Boost_FOUND)
    # If config mode fails, try using FindBoost as fallback
    # Set policy to suppress the warning about FindBoost being deprecated
    if(POLICY CMP0167)
        cmake_policy(SET CMP0167 OLD)
    endif()
    
    message(STATUS "Trying fallback method to find Boost...")
    find_package(Boost 1.87.0)
endif()

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost 1.87.0 or newer not found. Build cannot continue. Please install Boost 1.87.0 or newer.")
endif()

message(STATUS "Boost found at: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost version: ${Boost_VERSION}")
include_directories(${Boost_INCLUDE_DIRS})
add_compile_definitions(HAVE_BOOST)

# Find and enable OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found. Adding compiler flags: ${OpenMP_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


FetchContent_Declare(matplotplusplus
        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
        GIT_TAG origin/master)

FetchContent_GetProperties(matplotplusplus)
if(NOT matplotplusplus_POPULATED)
    FetchContent_Populate(matplotplusplus)
    add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# Include directory
include_directories(${CMAKE_SOURCE_DIR}/include)

# Collect all source and header files
file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src/*.cc)
file(GLOB_RECURSE HEADERS ${CMAKE_SOURCE_DIR}/include/*.h)

# Create executable or library
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Link OpenMP
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)


# Link Matplot++ (Even though it's header-only, this ensures compatibility)
target_link_libraries(${PROJECT_NAME} PUBLIC matplot)

