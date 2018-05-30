include_directories(../mco/include)
include_directories(../tclap)
include_directories(../cli_lib/include)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})
include_directories(${GMP_INCLUDE_PATH})

set(SOURCE_FILES
main.cpp

#### Core
basic/modules.h
basic/modules.cpp

#### Modules

## EST
modules/est_benson_module.h
modules/est_benson_module.cpp
modules/est_ok10_module.h
modules/est_ok10_module.cpp

## AP
modules/ap_bf_module.h
modules/ap_bf_module.cpp
modules/ap_benson_module.h
modules/ap_benson_module.cpp
modules/ap_ok10_module.h
modules/ap_ok10_module.cpp

## EP
modules/ep_benson_module.h
modules/ep_benson_module.cpp
modules/ep_bs_module.h
modules/ep_bs_module.cpp
modules/ep_ls_module.h
modules/ep_ls_module.cpp
modules/ep_bs_par_module.h
modules/ep_bs_par_module.cpp
modules/ep_martins_module.h
modules/ep_martins_module.cpp
modules/ep_lc_approx_module.h
modules/ep_lc_approx_module.cpp
modules/ep_tsaggouris_module.h
modules/ep_tsaggouris_module.cpp
modules/ep_ideal_point.h
modules/ep_ideal_point.cpp
modules/ep_lex_dijkstra.h
modules/ep_lex_dijkstra.cpp
modules/ui_verifier.h
modules/ui_verifier.cpp
        basic/FrontWriter.h)

add_definitions(-DGMPRATIONAL)

add_executable(mco_cli ${SOURCE_FILES})

target_link_libraries(mco_cli mco)
target_link_libraries(mco_cli debug ${OGDF-DBG} optimized ${OGDF})
target_link_libraries(mco_cli debug ${COIN-DBG} optimized ${COIN})
target_link_libraries(mco_cli cdd)
target_link_libraries(mco_cli ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(mco_cli ${GMP})
target_link_libraries(mco_cli ${GMPXX})
target_link_libraries(mco_cli ${PPL})