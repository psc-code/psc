
#pragma once

#include <iterator>

template<typename T>
struct Range
{
  struct iterator : std::iterator<std::input_iterator_tag, T>
  {
    explicit iterator(T val) : val_(val) {}
    
    T operator*() const { return val_; }
    T const* operator->() const { return &val_; }

    iterator& operator++() { ++val_; return *this; }
    iterator operator++(int) { auto retval = *this; ++*this; return retval; }

    bool operator==(const iterator& other) const { return val_ == other.val_; }
    bool operator!=(const iterator& other) const { return !(*this == other);  }

  private:
    T val_;
  };

  Range(T begin, T end)
    : begin_(begin), end_(end)
  {}

  iterator begin() const { return begin_; }
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
    explicit iterator(T val, T stride) : val_(val), stride_(stride) {}
    
    T operator*() const { return val_; }
    T const* operator->() const { return &val_; }

    iterator& operator++() { val_ += stride_; return *this; }
    iterator operator++(int) { auto retval = *this; ++*this; return retval; }

    bool operator!=(const iterator& other) const { return val_ < other.val_; }

  private:
    T val_;
    T stride_;
  };

  RangeStrided(T begin, T end, T stride)
    : begin_(begin, stride), end_(end, stride)
  {}

  iterator begin() const { return begin_; }
  iterator end() const { return end_; }
  
private:
  iterator begin_;
  iterator end_;
};

template<typename T>
Range<T> range(T begin, T end)
{
  return {begin, end};
}

template<typename T>
Range<T> range(T length)
{
  return {0, length};
}

template<typename T>
RangeStrided<T> range(T begin, T end, T stride)
{
  return {begin, end, stride};
}

