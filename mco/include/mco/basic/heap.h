/*
 * This file is part of the paths library that is used to solve path problems in graphs.
 * Copyright (C) 2015  Denis Kurz <denis.kurz@tu-dortmund.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef PATHS_HEAP_H_
#define PATHS_HEAP_H_

#include <algorithm>
#include <functional>
#include <iterator>
#include <type_traits>

namespace mco
{

template<unsigned int D, class RandomIterator>
RandomIterator heap_parent(RandomIterator first, RandomIterator child)
{
    static_assert(std::is_same<std::random_access_iterator_tag,
                  typename std::iterator_traits<RandomIterator>::iterator_category>::value,
                 "Function heap_parent() only accepts random access iterators.");
    return first + (std::distance(first, child) - 1) / D;
}

template<unsigned int D, class RandomIterator>
RandomIterator first_heap_child(RandomIterator first, RandomIterator parent)
{
    static_assert(std::is_same<std::random_access_iterator_tag,
                  typename std::iterator_traits<RandomIterator>::iterator_category>::value,
                 "Function first_heap_child() only accepts random access iterators.");
    auto currentIndex = std::distance(first, parent);
    return first + (currentIndex * D) + 1;
}

template<unsigned int D, class RandomIterator>
std::size_t heap_child_count(RandomIterator first, RandomIterator it, RandomIterator end)
{
    using diff_type = typename RandomIterator::difference_type;
    auto zero = static_cast<diff_type>(0);
    auto d = static_cast<diff_type>(D);
    auto distToEnd = std::distance(first_heap_child<D>(first, it), end);
    return std::max(zero, std::min(distToEnd, d));
}

template<unsigned int D, class RandomIterator, class Compare>
void sift_down(RandomIterator first, RandomIterator it, RandomIterator last, Compare comp)
{
    static_assert(std::is_same<std::random_access_iterator_tag,
                  typename std::iterator_traits<RandomIterator>::iterator_category>::value,
                 "Function sift_down() only accepts random access iterators.");
    auto firstChild = first_heap_child<D>(first, it);
    auto lastChild = firstChild + D;
    while (lastChild <= last) {
        auto maxChild = std::max_element(firstChild, lastChild, comp);
        if (comp(*maxChild, *it)) return;
        std::iter_swap(it, maxChild);

        it = maxChild;
        firstChild = first_heap_child<D>(first, it);
        lastChild = firstChild + D;
    }
    if (firstChild < last) {
        auto maxChild = std::max_element(firstChild, last, comp);
        if (comp(*it, *maxChild)) {
            std::iter_swap(it, maxChild);
        }
    }
}

template<unsigned int D, class RandomIterator, class Compare>
bool is_ary_heap(RandomIterator first, RandomIterator last, Compare comp)
{
    bool result = true;
    auto child = first;
    ++child;
    while (child < last) {
        if (comp(*heap_parent<D>(first, child), *child)) {
            result = false;
            break;
        }
        ++child;
    }
    return result;
}
template<unsigned int D, class RandomIterator>
bool is_ary_heap(RandomIterator first, RandomIterator last)
{
    using Compare = std::less<typename std::iterator_traits<RandomIterator>::value_type>;
    Compare comp;
    return is_ary_heap<D, RandomIterator, Compare>(first, last, comp);
}

template<unsigned int D, class RandomIterator, class Compare>
void make_ary_heap(RandomIterator first, RandomIterator last, Compare comp)
{
    static_assert(std::is_same<std::random_access_iterator_tag,
                  typename std::iterator_traits<RandomIterator>::iterator_category>::value,
                 "Function make_ary_heap() only accepts random access iterators.");
    for (auto i = std::distance(first, last) / 2; i >= 0; --i) {
        sift_down<D>(first, first + i, last, comp);
    }
}

template<unsigned int D, class RandomIterator>
void make_ary_heap(RandomIterator first, RandomIterator last)
{
    using Compare = std::less<typename std::iterator_traits<RandomIterator>::value_type>;
    Compare comp;
    make_ary_heap<D, RandomIterator, Compare>(first, last, comp);
}

template<unsigned int D, class RandomIterator, class Compare>
void push_ary_heap(RandomIterator first, RandomIterator last, Compare comp)
{
    static_assert(std::is_convertible<
                    typename std::iterator_traits<RandomIterator>::iterator_category,
                    std::bidirectional_iterator_tag>::value,
                 "Function push_ary_heap() only accepts bidirectional iterators.");
    auto cur = last;
    --cur;
    auto parent = heap_parent<D>(first, cur);
    while (first < cur && comp(*parent, *cur)) {
        std::iter_swap(parent, cur);
        cur = parent;
        parent = heap_parent<D>(first, cur);
    }
}

template<unsigned int D, class RandomIterator>
void push_ary_heap(RandomIterator first, RandomIterator last)
{
    using Compare = std::less<typename std::iterator_traits<RandomIterator>::value_type>;
    Compare comp;
    push_ary_heap<D, RandomIterator, Compare>(first, last, comp);
}

template<unsigned int D, class RandomIterator, class Compare>
void pop_ary_heap(RandomIterator first, RandomIterator last, Compare comp)
{
    static_assert(std::is_convertible<
                    typename std::iterator_traits<RandomIterator>::iterator_category,
                    std::bidirectional_iterator_tag>::value,
                 "Function pop_ary_heap() only accepts bidirectional iterators.");
    --last;
    std::iter_swap(first, last);
    sift_down<D>(first, first, last, comp);
}

template<unsigned int D, class RandomIterator>
void pop_ary_heap(RandomIterator first, RandomIterator last)
{
    using Compare = std::less<typename std::iterator_traits<RandomIterator>::value_type>;
    Compare comp;
    pop_ary_heap<D, RandomIterator, Compare>(first, last, comp);
}

template<unsigned int D, class KeyType, class ValueType, class KeyCompare = std::less<KeyType>,
         class Container = std::vector<std::pair<KeyType, ValueType>>>
class MaxAryHeap
{

    static_assert(std::is_same<std::random_access_iterator_tag,
                  typename std::iterator_traits<typename Container::iterator>::iterator_category>::value,
                  "Class MaxAryHeap only accepts random access iterators.");
public:

    using Pair = std::pair<KeyType, ValueType>;

    MaxAryHeap(KeyCompare const &compare = KeyCompare())
    : compare_(compare)
    {
        pairCompare_ = [this](Pair const &p1, Pair const &p2) {
            return compare_(p1.first, p2.first);
        };
    }

    MaxAryHeap(MaxAryHeap const &other)
    : container_(other.container_), compare_(other.compare_), pairCompare_(other.pairCompare_)
    {
    }

    MaxAryHeap(MaxAryHeap &&other) noexcept
    : MaxAryHeap()
    {
        swap(*this, other);
    }

    MaxAryHeap &operator=(MaxAryHeap other)
    {
        swap(*this, other);
        return *this;
    }

    template<unsigned int Dprime, class KT, class VT, class KC>
    friend void swap(MaxAryHeap<Dprime, KT, VT, KC> &h1,
                     MaxAryHeap<Dprime, KT, VT, KC> &h2);

    void push(KeyType const &key, ValueType const &value)
    {
        auto pair = std::make_pair(key, value);
        container_.push_back(std::move(pair));
        push_ary_heap<D>(container_.begin(), container_.end(), pairCompare_);
    }

    ValueType pop()
    {
        auto result = container_[0].second;
        pop_ary_heap<D>(container_.begin(), container_.end(), pairCompare_);
        container_.pop_back();
        return result;
    }

    ValueType const &top()
    {
        return container_[0].second;
    }

    KeyType const &top_key()
    {
        return container_[0].first;
    }

    bool empty()
    {
        return container_.empty();
    }

    const Container& data()
    {
        return container_;
    }

private:

    Container container_;

    KeyCompare compare_;

    std::function<bool(Pair const &, Pair const &)> pairCompare_;

};

template<unsigned int D, class KeyType, class ValueType, class KeyCompare = std::less<KeyType>,
         class Container = std::vector<std::pair<KeyType, ValueType>>>
inline void swap(MaxAryHeap<D, KeyType, ValueType, KeyCompare, Container> &h1,
                 MaxAryHeap<D, KeyType, ValueType, KeyCompare, Container> &h2)
{
    using std::swap;
    swap(h1.container_, h2.container_);
    swap(h1.compare_, h2.compare_);
    swap(h1.pairCompare_, h2.pairCompare_);
}

} // namespace paths

#endif /* PATHS_HEAP_H_ */
