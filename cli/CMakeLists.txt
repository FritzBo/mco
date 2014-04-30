include_directories(../include)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})

set(SOURCE_FILES
main.cpp
)

#### Boost program options
find_library(BOOST_PO_LIBRARY
    NAMES   libboost_program_options-mt.a
    PATHS   /opt/local/
            /usr/local/
            /usr/
)

find_path(BOOST_PO_INCLUDE
    NAMES   boost/program_options.hpp
    PATHS   /opt/local/
            /usr/local/
            /usr/
)

include_directories(${BOOST_PO_INCLUDE})

add_executable(mco_cli ${SOURCE_FILES})

target_link_libraries(mco_cli mco)
target_link_libraries(mco_cli debug ${OGDF-DBG} optimized ${OGDF})
target_link_libraries(mco_cli debug ${COIN-DBG} optimized ${COIN})
target_link_libraries(mco_cli ${GUROBI1})
target_link_libraries(mco_cli ${GUROBI2})
target_link_libraries(mco_cli ${BOOST_PO_LIBRARY})
target_link_libraries(mco_cli ${CDD})
target_link_libraries(mco_cli pthread)