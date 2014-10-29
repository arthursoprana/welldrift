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

#ifndef ITL_NUMBER_TRAITS_H
#define ITL_NUMBER_TRAITS_H

#include <complex>

namespace itl {

  template <class T>
  struct number_traits {
    typedef T magnitude_type;
  };
  
#if defined _MSC_VER && !defined  __MWERKS__ && ! defined __ICL
  
  template <>
  struct number_traits<std::complex<float> > {
    typedef float magnitude_type;
  };
  template <>
  struct number_traits<std::complex<double> > {
    typedef double magnitude_type;
  };

#else
  template <class T>
  struct number_traits<std::complex<T> > {
    typedef T magnitude_type;
  };
#endif
  
}
#endif  //ITL_NUMBER_TRAITS_H
