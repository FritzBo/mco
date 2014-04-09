include_directories(../include)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})

set(SOURCE_FILES
main.cpp
)

add_executable(mco_cli ${SOURCE_FILES})

target_link_libraries(mco_cli mco)
target_link_libraries(mco_cli debug ${OGDF-DBG} optimized ${OGDF})
target_link_libraries(mco_cli debug ${COIN-DBG} optimized ${COIN})
target_link_libraries(mco_cli ${GUROBI})
target_link_libraries(mco_cli pthread)