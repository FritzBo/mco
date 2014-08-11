#pragma once
/*
 * product_iterator.h
 *
 *  Created on: 08.04.2013
 *      Author: fritz
 */

#ifndef PRODUCT_H_
#define PRODUCT_H_

#include <utility>
#include <iterator>
#include <vector>

#include <mco/point.h>

namespace mco {
    
template<typename T>
class Range {
    using value_type = T;
    using value_type_reference = T&;
    using const_value_type_reference = const T&;
    
public:
    Range(const_value_type_reference first,
          const_value_type_reference last)
    
    :   first_(first),
        last_(last) { }
    
    class Iterator {
    public:
        Iterator& operator++() {
            if(!end_ && value_type != range.last_) {
                value_type++;
            } else {
                end_ = true;
            }
            return *this;
        }
        const_value_type_reference operator*() const {
            return current_;
        }
        
        bool operator==(Iterator& iter) const {
            return this->current_ == iter.current_ &&
                this->end_ == iter.end_;
        }
        
        bool operator!=(Iterator& iter) const {
            return !this>operator==(iter);
        }
        
    private:
        Iterator(Range& range)
        :   range_(range),
            current_(range.first_),
            end_(false) { }
        
        Iterator(Range& range, const_value_type_reference current, bool end)
        :   range_(range),
            current_(current),
            end_(end) { }
        
        Range& range_;
        const_value_type_reference current_;
        bool end_;
    };
    
    Iterator begin() {
        return Iterator(*this);
    }
    
    Iterator end() {
        return Iterator(*this, last_, true);
    }
    
    pair<Iterator, Iterator> iterator_pair() {
        return make_pair(begin(), end());
    }
    
private:
    const_value_type_reference first_;
    const_value_type_reference last_;
};

template<typename ConstIterator>
class Product {
    
    using value_type = ConstIterator::value_type;

public:
	Product(std::vector<std::pair<ConstIterator, ConstIterator>> ranges)
    :   ranges_(ranges),
        begin_(begin),
        end_(false) {
            
		number_ranges_ = ranges.size();
		current_.reserve(number_ranges_);
		for(unsigned int i = 0; i < number_ranges_; ++i)
			current_[i] = ranges_[i].first;
	}
    
    class Iterator {
    public:
        Iterator& operator++() {
            
        }
        
        
        
        
        
    };
    
    Iterator begin() {
        
    }
    
    Iterator end() {
        
    }

	void operator()() {

		if(end_ == true)
			return;

		OutputIterator output(begin_);

		for(auto iterator : current_) {
			*output = *iterator;
			output++;
		}

		for(unsigned int i = 0; i < number_ranges_; ++i) {
			if(current_[i] != ranges_[i].second) {
				current_[i]++;
				break;
			} else {
				current_[i] = ranges_[i].first;
				if(i == number_ranges_ - 1)
					end_ = true;
			}
		}
	}

private:
    
    std::vector<std::pair<InputIterator, InputIterator>> ranges_;
	OutputIterator begin_;
	std::vector<InputIterator> current_;
	unsigned int number_ranges_;
	bool end_;

};

template<typename InputIterator, typename OutputIterator>
Product<InputIterator, OutputIterator> generate_product(std::vector<std::pair<InputIterator, InputIterator>> ranges, OutputIterator begin) {
	return Product<InputIterator, OutputIterator>(ranges, begin);
}

}

#endif /* PRODUCT_H_ */
