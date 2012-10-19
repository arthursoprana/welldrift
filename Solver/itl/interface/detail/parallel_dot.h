#ifndef PARALLEL_DOT_H
#define PARALLEL_DOT_H

#include <mpi.h>
#include "Manager.h"

template <class VectorX, class VectorY>
inline double parallel_dot(const VectorX& x, const VectorY& y) {

  double local, g = 0;
  local = mtl::dot(x, y);
  
  MPI_Allreduce(&local, &g, 1, MPI_DOUBLE, MPI_SUM, Manager::comm());
  return g;
}

template <class VectorX, class VectorY>
inline double parallel_dot_conj(const VectorX& x, const VectorY& y) {

  double local, g = 0;
  local = mtl::dot_conj(x, y);
  
  MPI_Allreduce(&local, &g, 1, MPI_DOUBLE, MPI_SUM, Manager::comm());
  return g;
}


#endif // PARALLEL_DOT_H
