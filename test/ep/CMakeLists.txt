include_directories(../../include)
include_directories(../include)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})
include_directories(${GTEST_INCLUDE_PATH})
include_directories(${CDD_INCLUDE_PATH})

set(SOURCE_FILES
../include/instance_fixture.h

ep_benson_dual_test.cpp
lagrange_bisect_test.cpp
)

add_executable(ep_test ${SOURCE_FILES})

target_link_libraries(ep_test mco)
target_link_libraries(ep_test debug ${OGDF-DBG} optimized ${OGDF})
target_link_libraries(ep_test debug ${COIN-DBG} optimized ${COIN})
target_link_libraries(ep_test pthread)
target_link_libraries(ep_test ${CDD})
target_link_libraries(ep_test ${GTEST} ${GTEST_MAIN})

add_test(Efficient_Paths ep_test)