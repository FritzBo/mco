//
//  main.cpp
//  cli
//
//  Created by Fritz Bökler on 25.04.2016.
//
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <cctype>

using namespace std;

vector<vector<double>> read_points(string front_file_name)
{
    ifstream front_file(front_file_name);

    if(!front_file.good())
    {
        throw exception();
    }

    vector<vector<double>> point_container;

    front_file >> noskipws;

    double entry;
    while(front_file.peek() != EOF)
    {
        vector<double> point;
        while(front_file.peek() != '\n' &&
              front_file.peek() != EOF)
        {
            if(std::isdigit(front_file.peek()))
            {
                front_file >> entry;
                point.push_back(entry);
            }
            else
            {
                front_file.ignore();
            }
        }
        front_file.ignore();
        point_container.push_back(point);
    }

    if(point_container[point_container.size() - 1].size() == 0)
    {
        point_container.pop_back();
    }

    return point_container;
}

double compute_ratio(string front_file_name,
                     string approx_file_name)
{
    vector<vector<double>> front = read_points(front_file_name);
    vector<vector<double>> approx = read_points(approx_file_name);

    if(front.size() == 0 || approx.size() == 0)
    {
        return numeric_limits<double>::infinity();
    }

    unsigned components = front[0].size();

    double ratio = 1;

    for(auto& front_point : front)
    {
        double point_ratio = numeric_limits<double>::infinity();

        for(auto& approx_point : approx)
        {
            double comp_ratio = 1;
            for(unsigned i = 0; i < components; ++i)
            {
                if(front_point[i] > 0)
                {
                    comp_ratio = max(comp_ratio, approx_point[i] / front_point[i]);
                }
            }
            point_ratio = min(point_ratio, comp_ratio);
        }
        ratio = max(ratio, point_ratio);
    }
    return ratio;
}

int main(int argc, char** argv)
{
    // todo: remove
#ifndef NDEBUG
    this_thread::sleep_for(std::chrono::seconds(1));
#endif

    if(argc != 3)
    {
        cout << "Syntax: " << argv[0] << "<front file> <approx file>" << endl;
        return 0;
    }

    string front_file_name(argv[1]);
    string approx_file_name(argv[2]);

    try {
        cout << compute_ratio(front_file_name, approx_file_name) << endl;
    } catch(exception& e) {
        cout << "Unable to read file(s)" << endl;
    }

    return 0;
}