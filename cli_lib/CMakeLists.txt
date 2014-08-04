include_directories(include)
include_directories(../include)
include_directories(../tclap)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})

set(SOURCE_FILES
# SOURCES

#### Dummy
src/dummy.cpp

## HEADER

#### Output
include/mco/cli/output_formatter.h
)

add_library(cli_lib STATIC ${SOURCE_FILES})