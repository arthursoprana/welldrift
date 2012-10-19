#ifndef NTL_PARALLEL_MATRIX_H
#define NTL_PARALLEL_MATRIX_H

/*
  There are several ways to have parallel matrix-vector
  multiplication.
  
  1) the naive one is to send/recv all x elements and do 
     computation.
  
  2) other alternative is to figure out the exact access pattern.
     setup persistent communication.
     then, either each time matrix-vector multiplication happens,

     1> start communication and wait until all data available 
        and do cumputation. Or

     2> start up communication, 
        do diagonal_block * local_x 
        wait all data avaiable
	do all off-diagonal block * other_x, or

     3> start up communication, 
        do diagonal_block * local_x 
        wait any data avaiable
	do the off-diagonal block * remote_x_available
*/


#include <mpi.h>

#include "itl/interface/detail/Manager.h"
#include "itl/interface/detail/pointer_vector.h"
#include "itl/interface/detail/mtl_amend.h"

/* This version of sparse_matrix only send/recv sparse blocks of
   vector x during mat-vec operation A*x.
 */


#include "itl/interface/detail/pattern_finder.h"

extern "C" {
  void amux(int* n, double* x, double* y, double* a, int* ja, int* ia) {
    for(int i=0; i<*n; ++i) {
      double t = 0;
      int s = ia[i], e = ia[i+1];
      for(int j=s; j<e; ++j)
	t += a[j] * x[ ja[j] ];

      y[i] = t;
    }
  }

}


#define MV_ASKING_TAG  2000
#define MV_INDICES_TAG 2999
#define MV_DATA_TAG    3000

template <class Matrix>
inline int ncols(const Matrix& m)
{
  return m.ncols();
}
template <class Matrix>
inline int nrows(const Matrix& m)
{
  return m.nrows();
}

template <class Matrix>
class parallel_matrix {
public:
  typedef Matrix matrix_type;
  typedef typename Matrix::value_type value_type;
  
  parallel_matrix(const Matrix& AA)
    : A(AA), buffer(new value_type[ncols(A)]), wx(buffer, ncols(A)) {

    int processors = Manager::size();
    int rank = Manager::rank();
    int pos = wx.size() / processors; 

    request = new MPI_Request[processors*2];

    int* send_sizes = new int[processors];
    int* recv_sizes = new int[processors];
    for (int i=0; i<processors; i++)
      recv_sizes[i] = send_sizes[i] = 0;

    int* send_displacements = new int[wx.size()];
    int* recv_displacements = new int[wx.size()];
    for (int i=0; i<wx.size(); i++) 
      send_displacements[i] = recv_displacements[i] = 0;

    //post recv asking message
    int count = 0;
    for (int i=0; i<processors; i++) {
      if ( i != rank ) {
	MPI_Irecv( send_sizes+i, 1, MPI_INT, i, 
		   MV_ASKING_TAG, Manager::comm(), 
		   &request[count++]);
      }
    }

    //analyze matrix data and decide which blocks of x need fetch.
#if (defined ITL_USE_BLAS) || 1
    typename Matrix::const_iterator i=A.begin(), iend = A.end();
    typename Matrix::OneD::const_iterator j, jend;
    for (; i != iend; ++i) {
      j = (*i).begin(); jend = (*i).end();
      for (; j != jend; ++j) 
	recv_displacements[j.column()] = 1;
    }
#else
    pattern_finder pf(recv_displacements);
    mtl::mult(A, pf, pf);
#endif    
    //make index contiguous
    for (int i=0; i<processors; i++) {
      if ( i == rank ) continue;
      int b = i*pos;
      int k=0;
      int* rind = recv_displacements + b;

      for (int j=0; j<pos; j++)
	if ( rind[j] )
	  rind[k++] =  j; //local index

      recv_sizes[i] = k;
    }
    //leftover
    int b = (processors-1)*pos + recv_sizes[processors-1];
    int k = recv_sizes[processors-1];
    for (int i = processors*pos; i<wx.size(); i++)
      if ( recv_displacements[i] ) {
	recv_displacements[b++] = i - (processors-1)*pos; //global index
	k++;
      }
    recv_sizes[processors-1] = k;

    count          = processors-1;
    for (int i=0; i<processors; i++) {
      if ( i != rank ) {
	MPI_Isend( recv_sizes+i, 1, MPI_INT, i,
		   MV_ASKING_TAG, Manager::comm(), request + count++);
      }
    }

    //send recv_sizes, recv send_sizes
    MPI_Waitall(count, request, MPI_STATUSES_IGNORE);

    count = 0;
    for (int i=0; i<processors; i++) {
      if ( i != rank && send_sizes[i] ) {
	MPI_Irecv( send_displacements + i*pos, send_sizes[i], MPI_INT, i, 
		   MV_INDICES_TAG, Manager::comm(), 
		   &request[count++]);
      }
    }

    count          = processors-1;
    for (int i=0; i<processors; i++) {
      if ( i != rank && recv_sizes[i] ) {

	MPI_Isend( recv_displacements+i*pos , recv_sizes[i], MPI_INT, i,
		   MV_INDICES_TAG, Manager::comm(), request + count++);
      }
    }

    //complete send request, (blocks is the array to send 1 or 0)
    MPI_Waitall(count, request, MPI_STATUSES_IGNORE);

    int* blocks = new int [wx.size()-(processors-1)*pos];
    for (int i=0; i < wx.size()-(processors-1)*pos ; i++)
      blocks[i] = 1;

    MVtype = new MPI_Datatype [processors*2];
    count = 0;
    for (int i=0; i<processors; i++) 
      if ( i != rank && send_sizes[i] != 0 ) {
	MPI_Type_indexed(send_sizes[i], blocks, send_displacements+i*pos,
			 MPI_DOUBLE, & MVtype[count]);
	MPI_Type_commit(&MVtype[count]);
	count++;
      }

    for (int i=0; i<processors; i++) 
      if ( i != rank && recv_sizes[i] != 0 ) {
	MPI_Type_indexed(recv_sizes[i], blocks, recv_displacements+i*pos,
			 MPI_DOUBLE, & MVtype[count]);
	MPI_Type_commit(&MVtype[count]);
	count++;
      }

    count = 0;
    //persisitent send
    for (int i=0; i<processors; i++) 
      if ( i!=rank && send_sizes[i] !=0 ) {
	MPI_Send_init( buffer+rank*pos, 1, MVtype[count], i, 
		       MV_DATA_TAG, Manager::comm(),
		       request + count);
	count++;
      }

    //persistent recv
    for (int i=0; i<processors; i++)
      if ( i!=rank && recv_sizes[i] !=0 ) {
	MPI_Recv_init( buffer+i*pos, 1, MVtype[count], i, 
		       MV_DATA_TAG, Manager::comm(), 
		       request + count);
	count++;
      }

    total_num_send_recv = count;

    delete [] blocks;
    delete [] recv_displacements;
    delete [] send_displacements;
    delete [] send_sizes;
    delete [] recv_sizes;
  }

  template <class Vector>
  inline void mult(Vector& b) const {
#if defined ITL_USE_BLAS
    int nrows = A.nrows();
    amux(&nrows, buffer, b.data(), 
	 const_cast<double*>(A.get_val()), 
	 const_cast<int*>(A.get_ind()), 
	 const_cast<int*>(A.get_ptr()));
#else
    mtl::mult(A, wx, b);
#endif
  }

  template <class Vector>
  void process_vector(const Vector& x) const {
    int processors = Manager::size();
    int rank = Manager::rank();
    int pos  = wx.size() / processors;
#if defined ITL_USE_BLAS    
    pointer_vector<value_type> local(const_cast<value_type*>(buffer+pos*rank),
				     x.size());
    itl::copy(x, local);
#else
    std::copy(x.begin(), x.end(), buffer + pos*rank);
#endif

    MPI_Startall(total_num_send_recv, request);

    //wait all sending/recving data 
    MPI_Waitall(total_num_send_recv, request, MPI_STATUSES_IGNORE);

  }

  inline Matrix local_matrix() const {
    return A;
  }

  inline ~parallel_matrix() {
    
    delete [] buffer;
    delete [] request;
    delete [] MVtype;
  }

  void cleanup() {
    for (int i=0; i<total_num_send_recv; i++) {
      //need free MPI_request here
      MPI_Request_free(request+i);
      //need free MPI_Type
      MPI_Type_free(MVtype+i);
    }

  }
protected:
  const Matrix& A;
  mutable value_type* buffer;
  pointer_vector<value_type> wx;
  
  MPI_Request* request;
  MPI_Datatype* MVtype;
  int total_num_send_recv;
};

#endif
