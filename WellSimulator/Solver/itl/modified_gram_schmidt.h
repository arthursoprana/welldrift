// -*- c++ -*-
//
//=======================================================================
// Copyright (C) 1997-2001
// Authors: Andrew Lumsdaine <lums@osl.iu.edu> 
//          Lie-Quan Lee     <llee@osl.iu.edu>
//
// This file is part of the Iterative Template Library
//
// You should have received a copy of the License Agreement for the
// Iterative Template Library along with the software;  see the
// file LICENSE.  
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.
//=======================================================================
//

#ifndef ITL_MODIFIED_GRAM_SCHMIDT_ORTHOGONALIZATION_H
#define ITL_MODIFIED_GRAM_SCHMIDT_ORTHOGONALIZATION_H

#include <vector>

#include <itl/itl.h>
#include <itl/itl_tags.h>

namespace itl {

  template <class Vec>
  class modified_gram_schmidt {
  public:
    typedef typename itl_traits<Vec>::size_type size_type; 

    template <class Size>
    modified_gram_schmidt(int restart, const Size& s) 
      : V(restart+1) {
      typedef typename itl_traits<Vec>::vector_category VecCat;
      do_initialize(VecCat(), s);
    }
    
    const Vec& operator[](size_type i) const {
      return V[i];
    }

    Vec& operator[](size_type i) {
      return V[i];
    }
    
  protected:
    //such as mtl::dense1D
    template <class Size>
    void do_initialize(referencing_object_tag, const Size& s) {
      for (size_type i=0; i<V.size(); ++i)
	V[i] = Vec(s);
    }
    
    //such as blitz array
    template <class Size>
    void do_initialize(non_referencing_object_tag, const Size& s) {
      for (size_type i=0; i<V.size(); ++i)
	itl::resize(V[i], s);
    }
    
    std::vector<Vec> V;
  };
  


  template <class Vec, class VecHi, class Size>
  void orthogonalize(modified_gram_schmidt<Vec>& V,  const VecHi& _Hi, Size i)
 {
    VecHi& Hi = const_cast<VecHi&>(_Hi);
    
    for (Size k = 0; k <= i; k++) {
      Hi[k] = itl::dot_conj(V[i+1], V[k]);
      itl::add(itl::scaled(V[k], -Hi[k]), V[i+1]);
    }
  }
  
  template <class Vec, class VecS, class VecX, class Size>
  void combine(modified_gram_schmidt<Vec>& V, const VecS& s, VecX& x, Size i) {
    for (Size j = 0; j < i; ++j)
      itl::add(itl::scaled(V[j], s[j]), x);
  }
  
}

#endif
