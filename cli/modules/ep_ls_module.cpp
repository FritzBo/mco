//
//  ep_lc_approx_module.cpp
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#include "ep_ls_module.h"

#include <sstream>
#include <functional>

using std::map;
using std::string;
using std::list;
using std::vector;
using std::pair;
using std::function;
using std::make_pair;
using std::numeric_limits;
using std::cout;
using std::endl;
using std::chrono::steady_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::ArgException;
using TCLAP::ValueArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::SwitchArg;
using TCLAP::MultiArg;

#include <mco/benchmarks/temporary_graphs_parser.h>
#include <mco/basic/point.h>
#include <mco/ep/label_setting/label_setting.h>

using mco::TemporaryGraphParser;
using mco::Point;

using mco::ForwardStar;
using mco::node;
using mco::edge;
using mco::FSEdgeArray;
using mco::FSNodeArray;
using mco::ForwardStarFileReader;

using mco::MOSPLabelSetting;

void EpLsModule::perform(int argc, char** argv) {
    try {
        CmdLine cmd("Label Setting for the Multiobjective Shortest Path Problem.", ' ', "0.1");

        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");
        
        SwitchArg is_directed_arg("d", "directed", "Should the input be interpreted as a directed graph?", false);

        cmd.add(file_name_argument);
        cmd.add(is_directed_arg);

        cmd.parse(argc, argv);
        
        string file_name = file_name_argument.getValue();
        bool directed = is_directed_arg.getValue();

        ForwardStar graph;
        FSNodeArray<int> extern_ids(graph);
        FSEdgeArray<Point> raw_costs(graph);
        unsigned dimension;
        node source, target;
        
        ForwardStarFileReader file_reader;
        file_reader.read(file_name,
                         graph,
                         extern_ids,
                         raw_costs,
                         dimension,
                         source,
                         target,
                         directed);

        num_nodes_ = graph.numberOfNodes();
        num_edges_ = graph.numberOfEdges();
        num_objectives_ = dimension;

        auto cost_function = [&raw_costs] (edge e) -> Point* {
            return &raw_costs[e];
        };

        MOSPLabelSetting solver;

        steady_clock::time_point start = steady_clock::now();

        MOSPLabelSetting::InputDefinition def;

        def.graph = &graph;
        def.costs = cost_function;
        def.dimension = dimension;
        def.source = source;
        def.target = target;

        solver.Solve(def);

        steady_clock::time_point end = steady_clock::now();
        duration<double> computation_span
        = duration_cast<duration<double>>(end - start);
        solution_time_ = computation_span.count();

        for(auto& solution : solver.solutions())
        {
            list<mco::node> path;
            for(auto n : solution.first)
            {
                path.push_back(extern_ids[n]);
            }
            solutions_.push_back(make_pair(path, solution.second));
        }

        label_compares_ = 0;
        arc_pushes_ = 0;
        deleted_labels_ = 0;

    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const list<mco::node>, const Point>>& EpLsModule::
solutions() {
    return solutions_;
}

string EpLsModule::statistics() {
    std::stringstream sstream;
    sstream << num_nodes_;
    sstream << ", ";
    sstream << num_edges_;
    sstream << ", ";
    sstream << num_objectives_;
    sstream << ", ";
    sstream << label_compares_;
    sstream << ", ";
    sstream << deleted_labels_;
    sstream << ", ";
    sstream << arc_pushes_;
    sstream << ", ";
    sstream << solution_time_;

    return sstream.str();
}
