#ifndef NTL_PATTERN_FINDER_H
#define NTL_PATTERN_FINDER_H

class pattern_finder {
public:
  typedef int size_type;
  typedef int value_type;
  typedef const pattern_finder* const_iterator;
  typedef pattern_finder* iterator;
  inline pattern_finder(int* p_) : pattern(p_) {}

  inline int operator[](int i) const {
    pattern[i] = 1;
    return x;
  }
  const_iterator begin() const 
  { return this; }
  iterator begin() 
  { return this; }
protected:
  int x;
  mutable int* pattern;
};

struct pfinder_aux {
  inline int& operator[](int i) {
    return x;
  }
  int x;
};

#endif //NTL_PATTERN_FINDER_H
