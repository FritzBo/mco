include_directories(include)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})
include_directories(${GMP_INCLUDE_PATH})

set(SOURCE_FILES

#### Header ####

# Core
include/mco/basic/point.h
include/mco/basic/pareto_point_comparator.h
include/mco/basic/componentwise_point_comparator.h
include/mco/basic/equality_point_comparator.h
include/mco/basic/lex_point_comparator.h
include/mco/basic/abstract_solver.h
include/mco/basic/abstract_graph_instance.h
include/mco/basic/weight_function_adaptors.h
include/mco/basic/utility.h
include/mco/basic/cdd_files.h

# Assignment
include/mco/ap/basic/abstract_ap_solver.h
include/mco/ap/basic/ap_generator.h
include/mco/ap/basic/ap_instance.h
include/mco/ap/basic/lex_hungarian.h
include/mco/ap/benson_dual/ap_benson_dual_solver.h
include/mco/ap/benson_dual/ap_grb_benson_dual_solver.h
include/mco/ap/brute_force/ap_brute_force_solver.h
include/mco/ap/molp_solver/assignment_molp_model.h
include/mco/ap/ok10/ap_ok10_solver.h

# Benchmarks
include/mco/benchmarks/abstract_objective_vector_generator.h
include/mco/benchmarks/correlated_random.h
include/mco/benchmarks/graph_generators.h
include/mco/benchmarks/graph_parser.h
include/mco/benchmarks/kurz_parser.h
include/mco/benchmarks/mcap_parser.h
include/mco/benchmarks/pge10_parser3_ap.h
include/mco/benchmarks/random_objective_generator.h
include/mco/benchmarks/temporary_graphs_parser.h

# MO Shortest Path
include/mco/ep/basic/ep_instance.h
include/mco/ep/basic/abstract_ep_solver.h
include/mco/ep/basic/ep_lagrange_bisect.h
include/mco/ep/brum_shier/ep_solver_bs.h
include/mco/ep/brum_shier/ep_solver_bs_par.h
include/mco/ep/brum_shier/ep_weighted_bs.h
include/mco/ep/martins/label.h
include/mco/ep/martins/martins_smart.h
include/mco/ep/martins/martins.h
include/mco/ep/martins/weighted_martins.h
include/mco/ep/tsaggouris/ep_solver_tsaggouris_approx.h
include/mco/ep/warburton/ep_solver_warburton_approx.h
include/mco/ep/warburton/product.h
include/mco/ep/lc_approx/lc_approx.h
include/mco/ep/basic/binary_heap.h
include/mco/ep/basic/dijkstra.h
include/mco/ep/dual_benson/ep_dual_benson.h
include/mco/ep/dual_benson/constrained_ep_dual_benson.h
include/mco/ep/dual_benson/ep_benson_verifier.h
include/mco/ep/preprocessing/constrained_reach.h
include/mco/ep/basic/instance_scalarizer.h


# MO Spanning Tree
include/mco/est/basic/abstract_est_solver.h
include/mco/est/basic/kruskal_st_solver.h
include/mco/est/2tree/est_solver2_trees.h
include/mco/est/dual_benson/est_dual_benson_scalarizer.h

# Generics
include/mco/generic/lagrange_bisect.h

# Benson Dual
include/mco/generic/benson_dual/abstract_online_vertex_enumerator.h
include/mco/generic/benson_dual/dual_benson_scalarizer.h
include/mco/generic/benson_dual/ove_cdd.h
include/mco/generic/benson_dual/ove_cdd_gmp.h
include/mco/generic/benson_dual/ove_node_lists.h
include/mco/generic/benson_dual/ove_edge_lists.h
include/mco/generic/benson_dual/ove_fp_v2.h
include/mco/generic/benson_dual/upper_image_container.h

# Generics
include/mco/generic/ok10/ok10_scalarizer.h

# Geometry Tools
include/mco/geometric/projective_geometry_utilities.h

# MO Linear Programming
include/mco/molp/basic/molp_model.h
include/mco/molp/benson/abstract_vertex_container.h
include/mco/molp/benson/benson_exact.h
include/mco/molp/benson/benson.h
include/mco/molp/weightset/weightset_solver.h


#### Sources ####

# Assignment
src/ap/basic/ap_generator.cpp
src/ap/basic/lex_hungarian.cpp
src/ap/benson_dual/ap_benson_dual_solver.cpp
src/ap/benson_dual/ap_grb_benson_dual_solver.cpp
src/ap/brute_force/ap_brute_force_solver.cpp
src/ap/molp_solver/assignment_molp_model.cpp

# Benchmarks
src/benchmarks/correlated_random.cpp
src/benchmarks/graph_generators.cpp
src/benchmarks/graph_parser.cpp
src/benchmarks/kurz_parser.cpp
src/benchmarks/mcap_parser.cpp
src/benchmarks/pge10_parser3_ap.cpp
src/benchmarks/random_objective_generator.cpp
src/benchmarks/temporary_graphs_parser.cpp

# MO Shortest Path
src/ep/basic/ep_instance.cpp
src/ep/basic/ep_lagrange_bisect.cpp
src/ep/brum_shier/ep_solver_bs.cpp
src/ep/brum_shier/ep_solver_bs_par.cpp
src/ep/brum_shier/ep_weighted_bs.cpp
src/ep/martins/label.cpp
src/ep/martins/martins_smart.cpp
src/ep/martins/martins.cpp
src/ep/martins/weighted_martins.cpp
src/ep/tsaggouris/ep_solver_tsaggouris_approx.cpp
src/ep/warburton/ep_solver_warburton_approx.cpp
src/ep/lc_approx/lc_approx.cpp
src/ep/basic/dijkstra.cpp
src/ep/preprocessing/constrained_reach.cpp
src/ep/basic/instance_scalarizer.cpp

# MO Spanning Tree
src/est/basic/kruskal_st_solver.cpp
src/est/2tree/est_solver2_trees.cpp
src/est/dual_benson/est_dual_benson_scalarizer.cpp

# Generic
src/generic/lagrange_bisect.cpp

# Benson Dual
src/generic/benson_dual/ove_cdd.cpp
src/generic/benson_dual/ove_cdd_gmp.cpp
src/generic/benson_dual/ove_node_lists.cpp
src/generic/benson_dual/ove_edge_lists.cpp
src/generic/benson_dual/ove_fp_v2.cpp

# MO Linear Programming
src/molp/basic/molp_model.cpp
src/molp/benson/benson_exact.cpp
src/molp/benson/benson.cpp
src/molp/weightset/weightset_solver.cpp

# Doxygen
Doxyfile.in
)

add_definitions(-DGMPRATIONAL)

add_library(mco STATIC ${SOURCE_FILES})

## Doxygen

find_package(Doxygen)

if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile @ONLY)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)