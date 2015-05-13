cmake_minimum_required(VERSION 2.6)

project(mco)

set(VERSION_MAJOR   0   CACHE STRING "Project major version number.")
set(VERSION_MINOR   0   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   1   CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

enable_testing()

set(CMAKE_CXX_FLAGS
"--std=c++11 -pedantic -Wall -Wextra")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DOGDF_DEBUG")

## OGDF and COIN ##
find_library(OGDF
    NAMES libOGDF.a
    PATHS /usr/local/
	~/.local/
)

find_library(OGDF-DBG
    NAMES OGDF-dbg
    PATHS /usr/local/
)

find_path(OGDF_INCLUDE_PATH
    NAMES ogdf/basic/Graph.h
    PATHS /usr/local/
)

find_library(COIN
    NAMES COIN
    PATHS /usr/local/
)

find_library(COIN-DBG
    NAMES COIN-dbg
    PATHS /usr/local/
)

find_path(COIN_INCLUDE_PATH
    NAMES coin
    PATHS /usr/local/
)

## GUROBI ##

find_path(GUROBI_INCLUDE_PATH
    NAMES gurobi_c++.h
    PATHS /Library/*/mac64/include/
          /usr/local/
)

find_library(GUROBI1
    NAMES gurobi_c++
    PATHS /Library/*/mac64/lib/
          /usr/local/
)
          
find_library(GUROBI2
    NAMES gurobi56
    PATHS /Library/*/mac64/lib/
          /usr/local/
)

## PPL ##

find_path(PPL_INCLUDE_PATH
    NAMES   ppl.hh
    PATHS 	/usr/
            /usr/local/
            /opt/local/
)

find_library(PPL
    NAMES   libppl.a
    PATHS   /usr/
            /usr/local
            /opt/local
)

## GMP ##
find_path(GMP_INCLUDE_PATH
    NAMES   gmpxx.h
    PATHS 	/usr/
            /usr/local/
            /opt/local/
)

find_library(GMPXX
    NAMES   libgmpxx.a
    PATHS   /usr/
            /usr/local
            /opt/local
)

find_library(GMP
NAMES   libgmp.a
PATHS   /usr/
/usr/local
/opt/local
)

## Threads ##
find_package(Threads)

add_subdirectory(mco)
add_subdirectory(cdd)
add_subdirectory(test)
add_subdirectory(cli)
add_subdirectory(cli_lib)
add_subdirectory(cli_warburton)