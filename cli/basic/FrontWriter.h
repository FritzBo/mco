//
// Created by fritz on 30.05.18.
//

#ifndef MCO_FRONTWRITER_H
#define MCO_FRONTWRITER_H

#include <string>
#include <fstream>

#include <mco/basic/point.h>

class FrontWriter
{

public:
    FrontWriter() = default;

    template<typename OutputIterator>
    void WriteToFile(std::string, OutputIterator, OutputIterator);
};

template<typename OutputIterator>
void FrontWriter::
WriteToFile(std::string file_name, OutputIterator begin, OutputIterator end)
{
    std::ofstream file(file_name);

    OutputIterator iter = begin;

    while(iter != end)
    {
        mco::Point const & p = iter->second;
        for(unsigned i = 0; i < p.dimension() - 1; i++)
        {
            file << p[i] << ", ";
        }
        file << p[p.dimension() - 1] << std::endl;

        ++iter;
    }

    file.close();
}


#endif //MCO_FRONTWRITER_H
