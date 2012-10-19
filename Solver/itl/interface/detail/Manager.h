#ifndef NTL_MANAGER_H
#define NTL_MANAGER_H

#define MANAGER_MPI

#ifdef MANAGER_MPI
#include <mpi.h>
#endif

class Manager {
public:
  Manager() {}
  static void init(int& argc, char**& argv)
  {
#ifdef MANAGER_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &myComm);
    MPI_Comm_set_name(myComm, "ITL Comm");
    MPI_Comm_rank(myComm, &_rank);
    MPI_Comm_size(myComm, &_size);
#else
    _size = 1;
    _rank = 0;
#endif
  }

  static int size() { return _size; }
  static int rank() { return _rank; }
#ifdef MANAGER_MPI
  static MPI_Comm& comm() { return myComm; }
#endif

  static void finish() { 
#ifdef MANAGER_MPI
    MPI_Finalize();
#endif
  }

  static int _size;
  static int _rank; 
#ifdef MANAGER_MPI
  static MPI_Comm myComm;
#endif

};

int Manager::_size;
int Manager::_rank;

#ifdef MANAGER_MPI
MPI_Comm Manager::myComm;
#endif

#endif //
