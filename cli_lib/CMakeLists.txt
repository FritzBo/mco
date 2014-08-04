include_directories(../include)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})

set(SOURCE_FILES
#### Dummy
src/dummy.cpp
)

add_library(cli_lib STATIC ${SOURCE_FILES})