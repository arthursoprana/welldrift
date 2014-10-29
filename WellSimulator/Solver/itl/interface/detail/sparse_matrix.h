#ifndef ITL_INTERFACE_BLAS_INTERNAL_SPARSE_MATRIX_H
#define ITL_INTERFACE_BLAS_INTERNAL_SPARSE_MATRIX_H

namespace itl {

  //used in blas interface example
  //it is not responsible for managing memory
  template <class T>
  class sparse_matrix {
  public:
    //in amux and blas, it is int* instead of size_t!
    typedef int size_type;
    typedef T   value_type;
    sparse_matrix(size_type nr, size_type nc, size_type nnz__,
		  T* val_, size_type* ind_, size_type* ptr_) 
      : nrows_(nr), ncols_(nc), nnz_(nnz__),
      val(val_), ind(ind_), ptr(ptr_) {}
    const T* get_val() const { return val; }
    const size_type* get_ind() const { return ind; }
    const size_type* get_ptr() const { return ptr; }
    size_type nrows() const { return nrows_; }
    size_type ncols() const { return ncols_; }
    size_type nnz()   const { return nnz_; }
  private:
    size_type nrows_;
    size_type ncols_;
    size_type nnz_;
    T* val;
    size_type* ind;
    size_type* ptr;
  };


}


#endif //ITL_INTERFACE_BLAS_INTERNAL_SPARSE_MATRIX_H
