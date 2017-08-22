
cmake_minimum_required(VERSION 2.8)
project(ipopt_example)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(IPOPT_LIBRARY_DIRS /home/brian/Academics/Code_Projects/CoinIpopt/lib )
set(IPOPT_INCLUDE_DIRS /home/brian/Academics/Code_Projects/CoinIpopt/include/coin )
message("ipopt binary path is: " ${IPOPT_LIBRARY_DIRS} )
message("ipopt include path is: " ${IPOPT_INCLUDE_DIRS} )

include_directories(${IPOPT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
link_directories(${IPOPT_LIBRARY_DIRS})
# link_directories(${CMAKE_SOURCE_DIR})

add_library(opt_interface src/opt_interface.cpp)

add_executable(ipopt_example src/opt_main.cpp)
add_executable(unit_test_a src/unit_tests/unit_test_a.cpp)
add_executable(unit_test_b src/unit_tests/unit_test_b.cpp)

target_link_libraries(ipopt_example opt_interface)
target_link_libraries(ipopt_example ipopt)
