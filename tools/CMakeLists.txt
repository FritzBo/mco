include_directories(../mco/include)
include_directories(../tclap)
include_directories(${OGDF_INCLUDE_PATH})
include_directories(../cli_lib/include)

add_executable(appr-ratio src/appr_ratio.cpp)

add_executable(sp-tree src/sp-tree.cpp)
target_link_libraries(sp-tree mco)
target_link_libraries(sp-tree debug ${OGDF-DBG} optimized ${OGDF})