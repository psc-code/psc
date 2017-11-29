
#ifndef VPIC_LIST_BASE_H
#define VPIC_LIST_BASE_H

// ======================================================================
// VpicListBase

template<typename T>
class VpicListBase
{
public:
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;

  // ----------------------------------------------------------------------
  // class iterator
  
  class iterator
  {
    friend class VpicListBase<T>;

  protected:
    iterator(T* node) : node_(node) {}

  public:
    iterator() {}
    bool operator!=(const iterator& x) const { return node_ != x.node_; }
    reference operator*() const { return *node_; }
    pointer operator->() const { return node_; }
    iterator& operator++() { node_ = static_cast<T*>(node_->next); return *this; }

  private:
    T* node_;
  };

  // ----------------------------------------------------------------------
  // class const_iterator
  
  class const_iterator
  {
    friend class VpicListBase<T>;

  protected:
    const_iterator(const T* node) : node_(node) {}

  public:
    const_iterator() {}
    const_iterator(const iterator& x) : node_(x.node_) {}
    bool operator!=(const const_iterator& x) const { return node_ != x.node_; }
    const_reference operator*() const { return *node_; }
    const_pointer operator->() const { return node_; }
    const_iterator& operator++() { node_ = static_cast<T*>(node_->next); return *this; }

  private:
    const T* node_;
  };

  // ----------------------------------------------------------------------
  // VpicListBase itself

  VpicListBase() : head_(nullptr) {}

  bool empty() { return !head_; }
  
  iterator begin() { return iterator(head_); }
  iterator end()   { return iterator(nullptr); }

  const_iterator cbegin() const { return const_iterator(head_); }
  const_iterator cend()   const { return const_iterator(nullptr); }

protected:
  T* head_;
};

#endif

