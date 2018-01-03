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
#include <iterator>
#include <iomanip>

#include <cdd/setoper.h>
#include <cdd/cdd.h>

#include <mco/basic/point.h>

namespace mco {

template<typename OutputIterator>
class AbstractUpperImageContainer {
public:

    virtual OutputIterator cbeginExtremePoints() = 0;
    virtual OutputIterator cendExtremePoints() = 0;

    virtual OutputIterator cbeginInequalities() = 0;
    virtual OutputIterator cendInequalities() = 0;
};

class SimpleUpperImageContainer :
    public AbstractUpperImageContainer<std::list<Point*>::const_iterator> {

    using iterator = std::list<Point*>::const_iterator;

public:

    ~SimpleUpperImageContainer() {
        for(auto pp : extreme_points_list) {
            delete pp;
        }
        extreme_points_list.clear();
        for(auto pp : inequalities_list) {
            delete pp;
        }
        inequalities_list.clear();
    }

    std::list<Point*> extreme_points_list;
    std::list<Point*> inequalities_list;

    virtual iterator cbeginExtremePoints() {
        return extreme_points_list.cbegin();
    }
    virtual iterator cendExtremePoints() {
        return extreme_points_list.cend();
    }

    virtual iterator cbeginInequalities() {
        return inequalities_list.cbegin();
    }
    virtual iterator cendInequalities() {
        return inequalities_list.cend();
    }
};

template<typename OutputIterator>
class UpperImageWriter {
public:
    void write_image(std::string filename,
                     AbstractUpperImageContainer<OutputIterator>* container);

private:
    void write_geom_file(std::string filename,
                         OutputIterator cbegin,
                         OutputIterator cend);

};

template<typename T>
UpperImageWriter<T>
get_upper_image_writer(AbstractUpperImageContainer<T>* container) {

    return UpperImageWriter<T>();
}

template<typename OutputIterator>
void UpperImageWriter<OutputIterator>::
write_image(std::string filename,
            AbstractUpperImageContainer<OutputIterator>* container) {

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
        std::cout << "Could not open file: " << filename << std::endl;
        geom_file.close();
        return;
    }

    Point& test = **cbegin;
    unsigned dimension = test.dimension();

    unsigned distance = std::distance(cbegin, cend);

    geom_file << std::setprecision(100);

    geom_file << dimension << "\t" << distance << std::endl;

    auto it = cbegin;
    while(it != cend) {
        auto p = *it;
        for(unsigned i = 0; i < dimension; ++i) {
            geom_file << p->operator[](i);

            if(i < dimension - 1) {
                geom_file << "\t";
            }
        }
        
        ++it;

        if(it != cend) {
            geom_file << std::endl;
        }
    }

    geom_file.close();


}

class UpperImageReader {
public:
    void read_upper_image(std::string filename,
                          SimpleUpperImageContainer& container);

private:
    template<typename InputIterator>
    void read_geom_file(std::string filename,
                        InputIterator inserter);

};

inline void UpperImageReader::
read_upper_image(std::string filename,
                 SimpleUpperImageContainer& container) {

    read_geom_file(filename + ".bine", std::back_inserter(container.inequalities_list));

    read_geom_file(filename + ".bext", std::back_inserter(container.extreme_points_list));

}

template<typename InputIterator>
void UpperImageReader::
read_geom_file(std::string filename,
               InputIterator inserter) {

    std::fstream geom_file(filename, std::fstream::in);

    if(!geom_file.is_open()) {
        std::cout << "Could not open file: " << filename << std::endl;
        geom_file.close();
        return;
    }

    unsigned dimension, distance;
    geom_file >> dimension;
    geom_file >> distance;

    Point* p;
    for(unsigned j = 0; j < distance; ++j) {
        p = new Point(dimension);

        for(unsigned i = 0; i < dimension; ++i) {
            geom_file >> p->operator[](i);
        }

        inserter = p;

    }

    geom_file.close();

}

class UpperImageDDVerifier {

public:
    template<typename OutputIterator>
    void verify_double_description(AbstractUpperImageContainer<OutputIterator>& container);

private:
    template<typename OutputIterator>
    void init(AbstractUpperImageContainer<OutputIterator>& container);

    template<typename OutputIterator>
    void compare(AbstractUpperImageContainer<OutputIterator>& container);

    unsigned dimension_;
    dd_MatrixPtr h_representation_, v_representation_;
    dd_PolyhedraPtr polyhedron_;
};

template<typename OutputIterator>
void UpperImageDDVerifier::
verify_double_description(AbstractUpperImageContainer<OutputIterator>& container) {


    dd_set_global_constants();

    init(container);
    compare(container);

    dd_FreeMatrix(h_representation_);
    dd_FreeMatrix(v_representation_);
    dd_FreePolyhedra(polyhedron_);

    dd_free_global_constants();
    
}

template<typename OutputIterator>
void UpperImageDDVerifier::init(AbstractUpperImageContainer<OutputIterator>& container) {

    Point& test = **container.cbeginExtremePoints();
    dimension_ = test.dimension();

    unsigned number_ine = std::distance(
                                        container.cbeginInequalities(),
                                        container.cendInequalities());


    h_representation_ = dd_CreateMatrix(number_ine,
                                        dimension_ + 1);

    auto it = container.cbeginInequalities();
    unsigned line = 0;
    while(it != container.cendInequalities()) {

        long sum = 0;
        const long den = 1000000000000000;

        for(unsigned i = 0; i < dimension_ - 1; ++i) {
            long num = (*it)->operator[](i) * den;
            assert(num >= 0);

            dd_set_si2(h_representation_->matrix[line][i + 1],
                       num, den);

            sum += num;
        }

        long num = (*it)->operator[](dimension_ - 1) * den;
        dd_set_si2(h_representation_->matrix[line][0],
                   -num, den);

        long l_d = std::max((long) 0, den - sum);

        dd_set_si2(h_representation_->matrix[line][dimension_],
                 l_d, den);

        ++it;
        ++line;
    }
    
    h_representation_->representation = dd_Inequality;

    dd_ErrorType err;

    polyhedron_ = dd_DDMatrix2Poly(h_representation_, &err);

    if(err != dd_NoError) {
        dd_WriteErrorMessages(stdout, err);
    }

    v_representation_ = dd_CopyGenerators(polyhedron_);
}

template<typename OutputIterator>
void UpperImageDDVerifier::compare(AbstractUpperImageContainer<OutputIterator>& container) {

    unsigned num_rows = v_representation_->rowsize;

    std::cout << "Exact computation found " << num_rows << " extreme points." << std::endl;

//    dd_WriteMatrix(stdout, v_representation_);

    double max_pair_distance = - std::numeric_limits<double>::infinity();

    auto it = container.cbeginExtremePoints();
    while(it != container.cendExtremePoints()) {

        double min_distance = std::numeric_limits<double>::infinity();

        for(unsigned i = 0; i < v_representation_->rowsize; ++i) {

            if(dd_get_d(v_representation_->matrix[i][0]) == 1) {

                double distance = 0;
                for(unsigned j = 1; j < dimension_ + 1; ++j) {
                    distance += (dd_get_d(v_representation_->matrix[i][j]) - (*it)->operator[](j - 1)) * (dd_get_d(v_representation_->matrix[i][j]) - (*it)->operator[](j - 1));
                }

                distance = std::sqrt(distance);
                min_distance = std::min(min_distance, distance);
            }
        }

        ++it;

        max_pair_distance = std::max(max_pair_distance, min_distance);
    }

    std::cout << max_pair_distance << std::endl;

}

}

#endif
