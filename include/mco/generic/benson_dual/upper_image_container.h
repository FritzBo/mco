//
//  upper_image_container.h
//  mco
//
//  Created by Fritz Bökler on 10.11.14.
//
//

#ifndef mco_upper_image_container_h
#define mco_upper_image_container_h

#include <string>
#include <fstream>

#include <mco/basic/point.h>

namespace mco {

template<typename OutputIterator>
class UpperImageContainer {
public:

    virtual OutputIterator cbeginExtremePoints() = 0;
    virtual OutputIterator cendExtremePoints() = 0;

    virtual OutputIterator cbeginInequalities() = 0;
    virtual OutputIterator cendInequalities() = 0;
};

template<typename OutputIterator>
class UpperImageWriter {
public:
    void write_image(std::string filename,
                     UpperImageContainer<OutputIterator>* container);

private:
    void write_geom_file(std::string filename,
                   OutputIterator cbegin,
                   OutputIterator cend);

};

template<typename T>
UpperImageWriter<T>
get_upper_image_writer(UpperImageContainer<T>* container) {

    return UpperImageWriter<T>();
}

template<typename OutputIterator>
void UpperImageWriter<OutputIterator>::
write_image(std::string filename,
            UpperImageContainer<OutputIterator>* container) {

    write_geom_file(filename + ".bext",
                    container->cbeginExtremePoints(),
                    container->cendExtremePoints());

    write_geom_file(filename + ".bine",
                    container->cbeginInequalities(),
                    container->cendInequalities());

}

template<typename OutputIterator>
void UpperImageWriter<OutputIterator>::
write_geom_file(std::string filename,
                OutputIterator cbegin,
                OutputIterator cend) {

    std::fstream geom_file(filename, std::fstream::out);

    if(!geom_file.is_open()) {
        cout << "Could not open file: " << filename << endl;
        geom_file.close();
        return;
    }

    Point& test = **cbegin;

    geom_file << test.dimension() << std::endl;

    auto it = cbegin;
    while(it != cend) {
        auto p = *it;
        for(double val : *p) {
            geom_file << val << "\t";
        }
        geom_file << std::endl;
        
        it++;
    }

    geom_file.close();


}

}

#endif
