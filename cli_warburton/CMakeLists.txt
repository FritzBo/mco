include_directories(../mco/include)
include_directories(../tclap)
include_directories(../cli_lib/include)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})

set(SOURCE_FILES
main.cpp
)

add_executable(warburton ${SOURCE_FILES})

target_link_libraries(warburton mco)
target_link_libraries(warburton debug ${OGDF-DBG} optimized ${OGDF})
target_link_libraries(warburton debug ${COIN-DBG} optimized ${COIN})
target_link_libraries(warburton ${CMAKE_THREAD_LIBS_INIT})