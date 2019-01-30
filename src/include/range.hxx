
#pragma once

#include <iterator>

template<typename T>
struct Range
{
  struct iterator : std::iterator<std::input_iterator_tag, T>
  {
    __device__ __host__
    explicit iterator(T val) : val_(val) {}
    
    __device__ __host__
    T operator*() const { return val_; }
    __device__ __host__
    T const* operator->() const { return &val_; }

    __device__ __host__
    iterator& operator++() { ++val_; return *this; }
    __device__ __host__
    iterator operator++(int) { auto retval = *this; ++*this; return retval; }

    __device__ __host__
    bool operator==(const iterator& other) const { return val_ == other.val_; }
    __device__ __host__
    bool operator!=(const iterator& other) const { return !(*this == other);  }

  private:
    T val_;
  };

  __device__ __host__
  Range(T begin, T end)
    : begin_(begin), end_(end)
  {}

  __device__ __host__
  iterator begin() const { return begin_; }
  __device__ __host__
  iterator end() const { return end_; }
  
private:
  iterator begin_;
  iterator end_;
};

template<typename T>
struct RangeStrided
{
  struct iterator : std::iterator<std::input_iterator_tag, T>
  {
    __device__ __host__
    explicit iterator(T val, T stride) : val_(val), stride_(stride) {}
    
    __device__ __host__
    T operator*() const { return val_; }
    __device__ __host__
    T const* operator->() const { return &val_; }

    __device__ __host__
    iterator& operator++() { val_ += stride_; return *this; }
    __device__ __host__
    iterator operator++(int) { auto retval = *this; ++*this; return retval; }

    __device__ __host__
    bool operator!=(const iterator& other) const { return val_ < other.val_; }

  private:
    T val_;
    T stride_;
  };

  __device__ __host__
  RangeStrided(T begin, T end, T stride)
    : begin_(begin, stride), end_(end, stride)
  {}

  __device__ __host__
  iterator begin() const { return begin_; }
  __device__ __host__
  iterator end() const { return end_; }
  
private:
  iterator begin_;
  iterator end_;
};

template<typename T>
__device__ __host__
Range<T> range(T begin, T end)
{
  return {begin, end};
}

template<typename T>
__device__ __host__
Range<T> range(T length)
{
  return {0, length};
}

template<typename T>
__device__ __host__
RangeStrided<T> range(T begin, T end, T stride)
{
  return {begin, end, stride};
}

