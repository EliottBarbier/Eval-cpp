cmake_minimum_required(VERSION 3.0.0)
project(eval VERSION 0.1.0)

include(CTest)
enable_testing()

include_directories("./includes")

FILE(GLOB mainSources *.cpp)
FILE(GLOB cppSources cpp/*.cpp)

add_executable (eval
        ${mainSources}
        ${cppSources}
)

target_include_directories(eval PUBLIC
        "${PROJECT_SOURCE_DIR}/includes")
        

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
