//
//  benson_module.cpp
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#include "est_benson_module.h"

#include <set>

using std::set;
using std::string;
using std::list;
using std::pair;
using std::function;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::NodeArray;
using ogdf::node;
using ogdf::edge;

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::ArgException;
using TCLAP::ValueArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::SwitchArg;

#include <mco/benchmarks/kurz_parser.h>
#include <mco/basic/abstract_graph_instance.h>
#include <mco/est/dual_benson/est_dual_benson_scalarizer.h>
#include <mco/generic/benson_dual/ove_cdd.h>
#include <mco/generic/benson_dual/ove_cdd_gmp.h>
#include <mco/basic/point.h>

#include <mco/generic/benson_dual/upper_image_container.h>

using mco::Point;
using mco::KurzParser;
using mco::AbstractGraphInstance;
using mco::ESTDualBensonScalarizer;
using mco::OnlineVertexEnumeratorCDD;
using mco::OnlineVertexEnumeratorCddGmp;

using mco::UpperImageWriter;


void EstBensonModule::perform(int argc, char** argv) {
    try {
        CmdLine cmd("OK10 method to find the XND part of the Pareto-frontier of the efficient assignment problem.", ' ', "0.1");
        
        ValueArg<double> epsilon_argument("e", "epsilon", "Epsilon to be used in floating point calculations.", false, 1E-8, "epsilon");
        
        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");

        SwitchArg use_cdd_arg("", "cdd", "Using CDD Library for the vertex enumeration", false);

        SwitchArg use_cdd_gmp_arg("", "gmp", "Using CDD Library in GMP mode", false);

        SwitchArg use_gl_ove_arg("", "gl-ove", "Using the graphless online vertex enumerator", false);

        ValueArg<string> output_files_arg("o", "output", "Saves a discription of the upper image to .ine and .ext files.", false, "", "filname");


        cmd.add(epsilon_argument);
        cmd.add(file_name_argument);
        cmd.add(use_cdd_arg);
        cmd.add(use_cdd_gmp_arg);
        cmd.add(use_gl_ove_arg);
        cmd.add(output_files_arg);

        cmd.parse(argc, argv);
        
        string file_name = file_name_argument.getValue();
        string output_file_name = output_files_arg.getValue();
        double epsilon = epsilon_argument.getValue();
        bool use_cdd = use_cdd_arg.getValue();
        bool use_gmp = use_cdd_gmp_arg.getValue();
        bool use_gl_ove = use_gl_ove_arg.getValue();

        Graph graph;
        EdgeArray<Point *> edge_array(graph);

        KurzParser parser(file_name);

        parser.get_graph(graph,
                         edge_array);

        dimension = edge_array[graph.chooseEdge()]->dimension();
        nodes_ = graph.numberOfNodes();
        edges_ = graph.numberOfEdges();

        AbstractGraphInstance instance(graph, edge_array, dimension);

        if((instance.dimension() <= 4 && !use_cdd) || use_gl_ove) {
            ESTDualBensonScalarizer<> solver(epsilon);

            solver.Solve(instance);

            solutions_.insert(solutions_.begin(),
                              solver.solutions().cbegin(),
                              solver.solutions().cend());

            auto writer = get_upper_image_writer(&solver);
            writer.write_image(output_file_name, &solver);
        } else {
            if(use_gmp) {
                ESTDualBensonScalarizer<mco::OnlineVertexEnumeratorCddGmp> solver(epsilon);

                solver.Solve(instance);

                solutions_.insert(solutions_.begin(),
                                  solver.solutions().cbegin(),
                                  solver.solutions().cend());

                auto writer = get_upper_image_writer(&solver);
                writer.write_image(output_file_name, &solver);
            } else {
                ESTDualBensonScalarizer<mco::OnlineVertexEnumeratorCDD> solver(epsilon);

                solver.Solve(instance);

                solutions_.insert(solutions_.begin(),
                                  solver.solutions().cbegin(),
                                  solver.solutions().cend());

                auto writer = get_upper_image_writer(&solver);
                writer.write_image(output_file_name, &solver);
            }
        }

    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const list<edge>, const Point>>& EstBensonModule::solutions() {
    return solutions_;
}

string EstBensonModule::statistics() {
    string stats("");
    stats += to_string(dimension) + ", " + to_string(nodes_) + ", " + to_string(edges_);
    return stats;
}
