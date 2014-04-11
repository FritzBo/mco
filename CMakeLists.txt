cmake_minimum_required(VERSION 2.6)

project(mco)

enable_testing()

set(CMAKE_CXX_FLAGS
--std=c++11)

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

find_library(GUROBI
    NAMES gurobi_c++
    PATHS /Library/*/mac64/lib/
          /usr/local/
)

## CDD ##
find_path(CDD_INCLUDE_PATH
	NAMES cdd.h
	PATHS 	/usr/
		/usr/local/
		/opt/local/
)

find_library(CDD
	NAMES libcdd.a
	PATHS	/usr/
		/usr/local
		/opt/local
)

add_subdirectory(mco)
add_subdirectory(test)
add_subdirectory(cli)
