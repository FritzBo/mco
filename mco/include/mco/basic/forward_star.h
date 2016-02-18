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

using edge = unsigned;
using node = unsigned;

class ForwardStar;

class ForwardStarFileReader;

template<typename T>
class FSEdgeArray
{
public:

    FSEdgeArray() = delete;

    inline explicit FSEdgeArray(const ForwardStar& fs);

    inline FSEdgeArray(const ForwardStar& fs,
                       T object);

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
    const ForwardStar& ref_object_;

    friend ForwardStar;
    friend ForwardStarFileReader;
};

template<typename T>
class FSNodeArray
{
public:

    FSNodeArray() = delete;

    inline explicit FSNodeArray(const ForwardStar& fs);

    inline FSNodeArray(const ForwardStar& fs,
                       T object);

    inline T& operator[](node e)
    {
        return objects_[e];
    }

    inline const T& operator[](node e) const
    {
        return objects_[e];
    }

private:
    std::vector<T> objects_;
    const ForwardStar& ref_object_;

    friend ForwardStar;
    friend ForwardStarFileReader;
};

class ForwardStar
{
    class AdjacencyCollection;
    class NodeCollection;

    class NodeCollection
    {
        class NodeCollectionIterator
        {
        public:
            NodeCollectionIterator(unsigned begin)
            :   pos_(begin)
            {
            }

            inline node operator*()
            {
                return pos_;
            }

            inline NodeCollectionIterator& operator++()
            {
                ++pos_;
                return *this;
            }

            inline bool operator!=(NodeCollectionIterator& other)
            {
                return pos_ != other.pos_;
            }

            inline bool operator==(NodeCollectionIterator& other)
            {
                return pos_ == other.pos_;
            }

        private:
            unsigned pos_;
        };
    public:
        inline NodeCollection(ForwardStar& fs)
        :   fs_(fs)
        {
        }

        inline NodeCollectionIterator begin()
        {
            return NodeCollectionIterator(0);
        }

        inline NodeCollectionIterator end()
        {
            return NodeCollectionIterator(fs_.no_nodes);
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
        class AdjacencyCollectionIterator
        {
        public:
            AdjacencyCollectionIterator(unsigned begin)
            :   pos_(begin)
            {
            }

            inline edge operator*()
            {
                return pos_;
            }

            inline AdjacencyCollectionIterator& operator++()
            {
                ++pos_;
                return *this;
            }

            inline bool operator!=(AdjacencyCollectionIterator& other)
            {
                return pos_ != other.pos_;
            }

            inline bool operator==(AdjacencyCollectionIterator& other)
            {
                return pos_ == other.pos_;
            }

        private:
            unsigned pos_;
        };
    public:
        AdjacencyCollection(const ForwardStar& fs,
                            node n)
        :   fs_(fs),
            n_(n)
        {

        }

        inline AdjacencyCollectionIterator begin() const
        {
            return AdjacencyCollectionIterator(fs_.first_edge_[n_]);
        }

        inline AdjacencyCollectionIterator end() const
        {
            return AdjacencyCollectionIterator(fs_.first_edge_[n_ + 1]);
        }

    private:
        const ForwardStar& fs_;
        unsigned n_;
    };

    std::vector<unsigned> first_edge_;

    FSEdgeArray<node> heads_;
    FSEdgeArray<node> tails_;

    template<typename T>
    friend class FSEdgeArray;

    template<typename T>
    friend class FSNodeArray;

    friend class ForwardStarFileReader;
};

class ReverseStar
{

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
FSEdgeArray<T>::FSEdgeArray(const ForwardStar& fs)
:   objects_(fs.no_edges),
    ref_object_(fs)
{
}

template<typename T>
FSEdgeArray<T>::FSEdgeArray(const ForwardStar& fs,
                            T object)
:   ref_object_(fs),
    objects_(fs.no_edges,
             object)
{
}

template<typename T>
FSNodeArray<T>::FSNodeArray(const ForwardStar& fs)
:   objects_(fs.no_nodes),
    ref_object_(fs)
{
}

template<typename T>
FSNodeArray<T>::FSNodeArray(const ForwardStar& fs,
                            T object)
:   ref_object_(fs),
    objects_(fs.no_nodes,
             object)
{
}


}

#endif /* defined(__mco__forward_star__) */
