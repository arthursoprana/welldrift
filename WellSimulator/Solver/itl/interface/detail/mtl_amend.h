namespace mtl {

  template <class T>
    struct linalg_traits<pointer_vector<T> > {
      typedef oned_tag  dimension;
      typedef dense_tag sparsity;
      typedef T         value_type;
      typedef typename number_traits<value_type>::magnitude_type 
              magnitude_type;
    };

}
