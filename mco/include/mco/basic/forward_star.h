//
//  forward_star.h
//  mco
//
//  Created by Fritz Bökler on 09.01.16.
//
//

#ifndef __mco__forward_star__
#define __mco__forward_star__

#include <vector>
#include <string>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {

using edge = int;
using node = int;

constexpr int nulledge = -1;
constexpr int nullnode = -1;

class ForwardStar;

class ForwardStarFileReader;

class ReverseStarConstructor;

template<typename T>
class FSArray
{

public:

    FSArray() = delete;

    inline explicit FSArray(unsigned size)
    :   objects_(size)
    { }

    inline FSArray(unsigned size,
            T object)
    :   objects_(size, object)
    { }

    inline T& operator[](edge e)
    {
        return objects_[e];
    }

    inline const T& operator[](edge e) const
    {
        return objects_[e];
    }

private:
    std::vector<T> objects_;

    friend ForwardStar;
    friend ForwardStarFileReader;
    friend ReverseStarConstructor;
};

template<typename T>
class FSEdgeArray : public FSArray<T>
{
public:

    inline explicit FSEdgeArray(const ForwardStar& fs);

    inline FSEdgeArray(const ForwardStar& fs,
                       T object);

};

template<typename T>
class FSNodeArray : public FSArray<T>
{
public:

    FSNodeArray() = delete;

    inline explicit FSNodeArray(const ForwardStar& fs);

    inline FSNodeArray(const ForwardStar& fs,
                       T object);

};

class GraphObjectIterator
{
public:
    GraphObjectIterator(unsigned begin)
    :   pos_(begin)
    {
    }

    inline edge operator*()
    {
        return pos_;
    }

    inline GraphObjectIterator& operator++()
    {
        ++pos_;
        return *this;
    }

    inline bool operator!=(GraphObjectIterator& other)
    {
        return pos_ != other.pos_;
    }

    inline bool operator==(GraphObjectIterator& other)
    {
        return pos_ == other.pos_;
    }

private:
    unsigned pos_;
};


class ForwardStar
{
    class AdjacencyCollection;

    class NodeCollection
    {
    public:
        inline NodeCollection(ForwardStar& fs)
        :   fs_(fs)
        {
        }

        inline GraphObjectIterator begin() const
        {
            return GraphObjectIterator(0);
        }

        inline GraphObjectIterator end() const
        {
            return GraphObjectIterator(fs_.no_nodes);
        }
        
    private:
        ForwardStar& fs_;
    };

    unsigned no_nodes;
    unsigned no_edges;

public:

    ForwardStar()
    :   no_nodes(0),
        no_edges(0),
        nodes(*this),
        first_edge_(0),
        heads_(*this),
        tails_(*this)
    {

    }

    ForwardStar(unsigned no_nodes,
                unsigned no_edges)
    :   no_nodes(no_nodes),
        no_edges(no_edges),
        nodes(*this),
        first_edge_(no_nodes + 1),
        heads_(*this),
        tails_(*this)
    {
        first_edge_[no_nodes] = no_edges;
    }

    ForwardStar(ogdf::Graph& graph,
                bool directed);

    NodeCollection nodes;

    inline AdjacencyCollection adj_edges(node n) const
    {
        return AdjacencyCollection(*this, n);
    }

    inline node head(edge e) const
    {
        return heads_[e];
    }

    inline node tail(edge e) const
    {
        return tails_[e];
    }

    inline unsigned numberOfNodes() const
    {
        return first_edge_.size() - 1;
    }

    inline unsigned numberOfEdges() const
    {
        return no_edges;
    }


private:
    class AdjacencyCollection
    {
    public:
        AdjacencyCollection(const ForwardStar& fs,
                            node n)
        :   fs_(fs),
            n_(n)
        {

        }

        inline GraphObjectIterator begin() const
        {
            return GraphObjectIterator(fs_.first_edge_[n_]);
        }

        inline GraphObjectIterator end() const
        {
            return GraphObjectIterator(fs_.first_edge_[n_ + 1]);
        }

    private:
        const ForwardStar& fs_;
        unsigned n_;
    };

    void resize()
    {
        first_edge_.resize(no_nodes + 1);
        heads_.objects_.resize(no_edges);
        tails_.objects_.resize(no_edges);
    }

    std::vector<edge> first_edge_;

    FSEdgeArray<node> heads_;
    FSEdgeArray<node> tails_;

    template<typename T>
    friend class FSEdgeArray;

    template<typename T>
    friend class FSNodeArray;

    friend class ForwardStarFileReader;

    friend class ReverseStarConstructor;
};

using EdgeTrace = FSEdgeArray<edge>;

class ReverseStarConstructor
{
public:
    static void construct_reverse_star(const ForwardStar& fs,
                                       ForwardStar& rs,
                                       EdgeTrace& rs_fs_trace);

};

class ForwardStarFileReader
{
public:
    
    void read(std::string filename,
              ForwardStar& graph,
              FSNodeArray<int>& extern_node_ids,
              FSEdgeArray<Point>& weights,
              unsigned& dimension,
              node& source,
              node& target,
              bool directed);
};

template<typename T>
class ReverseEdgeArray {
public:
    inline explicit ReverseEdgeArray(const FSEdgeArray<T>& edge_array,
                                     const EdgeTrace& trace)
    :   edge_array_(edge_array),
        trace_(trace) { }

    inline T& operator[](edge e)
    {
        return edge_array_[trace_[e]];
    }

    inline const T& operator[](edge e) const
    {
        return edge_array_[trace_[e]];
    }

private:
    const FSEdgeArray<T>& edge_array_;
    const EdgeTrace& trace_;
};

template<typename T>
FSEdgeArray<T>::FSEdgeArray(const ForwardStar& fs)
:   FSArray<T>(fs.no_edges)
{
}

template<typename T>
FSEdgeArray<T>::FSEdgeArray(const ForwardStar& fs,
                            T object)
:   FSArray<T>(fs.no_edges,
               object)
{
}

template<typename T>
FSNodeArray<T>::FSNodeArray(const ForwardStar& fs)
:   FSArray<T>(fs.no_nodes)
{
}

template<typename T>
FSNodeArray<T>::FSNodeArray(const ForwardStar& fs,
                            T object)
:   FSArray<T>(fs.no_nodes,
               object)
{
}


}

#endif /* defined(__mco__forward_star__) */
