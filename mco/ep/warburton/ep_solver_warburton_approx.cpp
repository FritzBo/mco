/*
 * warburton.cpp
 *
 *  Created on: 02.04.2013
 *      Author: fritz
 */

#include <vector>
#include <limits>
#include <cmath>
#include <functional>
#include <list>

using std::vector;
using std::min;
using std::max;
using std::log;
using std::pow;
using std::ceil;
using std::floor;
using std::numeric_limits;
using std::function;
using std::list;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::NodeArray;
using ogdf::node;
using ogdf::edge;

#include <mco/ep/basic/dijkstra.h>
#include <mco/ep/warburton/ep_solver_warburton_approx.h>

namespace mco {

void EpSolverWarburtonApprox::Solve(ogdf::Graph& graph,
                                    std::function<Point*(edge)> cost_function,
                                    unsigned dimension,
                                    const ogdf::node source,
                                    const ogdf::node target,
                                    const Point &epsilon,
                                    list<Point> linear_bounds,
                                    unsigned int processes,
                                    double theta) {
    
	const unsigned int number_nodes = graph.numberOfNodes();

	vector<EdgeArray<double>> weight_functions(dimension - 1, EdgeArray<double>(graph, 0));
	NodeArray<double> distances(graph);
	NodeArray<edge> predecessor(graph);
	vector<double> min_e(dimension - 1, numeric_limits<double>::infinity());
	vector<double> max_e(dimension - 1, - numeric_limits<double>::infinity());
	vector<int> ub(dimension - 1);
	vector<int> lb(dimension - 1);
	vector<double> label_limits(dimension);
	edge e;
	double weight, d;
	int max_i, min_i;
    
    Dijkstra<double> sssp_solver;

	// Computing the bounds
	for(unsigned int k = 0; k < dimension - 1; ++k) {

		forall_edges(e, graph) {
			weight = (*cost_function(e))[k];
			min_e[k] = min(min_e[k], weight);
			max_e[k] = max(max_e[k], weight);
		}

		label_limits[k] = static_cast<int>(ceil((number_nodes - 1) * theta / epsilon[k]));

		max_i = static_cast<unsigned int>(ceil(log(min(max_e[k] * (number_nodes - 1), max_e[k] * (number_nodes - 1) / epsilon[k] ) ) / log(theta) ) + 1);
		min_i = static_cast<unsigned int>(log(max(1.0, min_e[k] * (number_nodes - 1) / epsilon[k] ) ) / log(theta) );

		for(int i = min_i; i < max_i; ++i) {
			d = static_cast<double>(pow(2, max_i - i));
			forall_edges(e, graph)
				weight_functions[k][e] = static_cast<int>(floor((*cost_function(e))[k] * (number_nodes - 1) / (epsilon[k] * d) ));

			// TODO: Error Message / Exception
            if(!sssp_solver.singleSourceShortestPaths(graph, weight_functions[k], source, predecessor, distances)) {
                
                cout << "Kein Weg gefunden!" << endl;
                return;
                
            }

			if(distances[target] > (number_nodes - 1) * theta / epsilon[k] || max_i - i == 1) {
                
				lb[k] = max_i - i + 1;
				ub[k] = static_cast<int>(ceil(log(min(max_e[k] * (number_nodes - 1), max_e[k] * (number_nodes - 1) / epsilon[k]) / log(theta)) + 1));
				break;
                
			}
		}
	}

	label_limits[dimension - 1] = numeric_limits<int>::infinity();

	cout << "bounds:" << endl;
	for(unsigned int k = 0; k < dimension - 1; ++k)
		cout << "k: " << k << ", lb: " << lb[k] << ", ub: " << ub[k] << ", limit: " << label_limits[k] << endl;

}

} // namespace mco
