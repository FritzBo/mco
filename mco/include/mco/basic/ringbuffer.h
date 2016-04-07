//
//  ringbuffer.h
//  mco
//
//  Created by Fritz Bökler on 05.04.16.
//
//

#ifndef mco_ringbuffer_h
#define mco_ringbuffer_h

namespace mco {

template<typename T>
class ring_buffer : private std::vector<T>
{
public:
    inline ring_buffer(unsigned capacity)
    :   std::vector<T>(capacity)
    {
    }

    inline void push_back(T object)
    {
        std::vector<T>::operator[](back_) = object;

        back_ = (back_ + 1) % std::vector<T>::size();
    }

    inline T front()
    {
        return std::vector<T>::operator[](front_);
    }

    inline void pop_front()
    {

#ifndef NDEBUG
        std::vector<T>::operator[](front_) = 0;
#endif

        front_ = (front_ + 1) % std::vector<T>::size();
    }

    inline bool empty()
    {
        return front_ == back_;
    }

private:
    unsigned front_ = 0;
    unsigned back_ = 0;
};

}

#endif
