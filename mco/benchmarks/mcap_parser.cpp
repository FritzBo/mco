/*
 * mcap_parser.cpp
 *
 *  Created on: 11.11.2013
 *      Author: fritz
 */

#include <mco/benchmarks/mcap_parser.h>

#include <set>
#include <list>

using std::make_shared;
using std::set;
using std::list;
using std::shared_ptr;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::node;
using ogdf::edge;

namespace mco {

shared_ptr<AssignmentInstance> MCAPParser::get_instance() {
	ifstream file(filename_);

	if(!file.good())
		    cerr << "Could not open file " << filename_ << endl;

	unsigned int num_agents, dim;

	file >> num_agents; //number of nodes
	file >> dim; 		//dimension of objective function;

	auto graph = make_shared<Graph>();
	auto edge_cost = make_shared<EdgeArray<Point *>>(*graph);
	auto agents = make_shared<set<node>>();
	list<node> jobs;

	for(unsigned int i = 0; i < num_agents; ++i) {
		node agent_node = graph->newNode();
		agents->insert(agent_node);
	}

	for(unsigned int i = 0; i < num_agents; ++i) {
		node job_node = graph->newNode();
		jobs.push_back(job_node);
		for(auto agent : *agents)
			graph->newEdge(agent, job_node);
	}

	double value;
	edge e;
	for(unsigned int i = 0; i < dim; ++i)
		for(auto agent : *agents)
			forall_adj_edges(e, agent){
				if(i == 0) {
					edge_cost->operator [](e) = new Point(0.0, dim);
				}
				file >> value;
				edge_cost->operator [](e)->operator [](i) = value;
			}

	auto instance = make_shared<AssignmentInstance>(graph, edge_cost, agents, dim);

	//close the file
	file.close();

	return instance;
}

} /* namespace mco */
