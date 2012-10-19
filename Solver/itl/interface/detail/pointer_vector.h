#ifndef POINTER_VECTOR_H
#define POINTER_VECTOR_H

//used in itl/interface/blas.h

//a wrapper for built-in array, it does not responsible for memory management
template <class T>
struct pointer_vector{
  //typedef mtl::oned_tag dimension;
  typedef T value_type;
  //typedef mtl::dense_tag sparsity;
  typedef T* iterator;
  typedef const T* const_iterator;
  
  inline pointer_vector() : buffer(NULL), n(0) {}
  inline pointer_vector(T* _b, int _n) : buffer(_b), n(_n) {}

  inline iterator begin() { return buffer; }
  inline iterator end() { return buffer+n; }

  inline const_iterator begin() const { return buffer; }
  inline const_iterator end() const { return buffer+n; }

  inline T& operator[](int i) { return buffer[i]; }
  inline const T& operator[](int i) const { return buffer[i]; }
  
  inline int size() const { return n; }
  inline void resize(int _n) {
    T* tmp = new T[_n];
    for (int i=0; i<n; ++i)
      tmp[i] = buffer[i];
    if (buffer)
      delete buffer;
    buffer = tmp;
    n = _n;
  }

  inline T* data() { return buffer; }
  inline const T* data() const { return buffer; }

  T* buffer;
  int n;
};

#endif //POINTER_VECTOR_H
