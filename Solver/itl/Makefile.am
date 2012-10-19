AUTOMAKE_OPTIONS = foreign

includedir	 = ${prefix}/include/itl

EXTRA_HEADERS	 =				\
		   Iteration_concept.h		\
		   Preconditioner_concept.h	\
		   Matrix_concept.h             \
		   givens_rotation.h            \
		   itl.h                        \
		   itl_utils.h			\
		   itl_tags.h			\
		   modified_gram_schmidt.h      \
		   matrix_free_operator.h       \
		   number_traits.h


EXTRA_DIST  =	$(EXTRA_HEADERS) krylov interface  preconditioner

make-install-dirs:
	-if test '!' -d $(includedir); then mkdir -p $(includedir); fi

install-data-hook: make-install-dirs
	-@ echo Installing $(includedir) ; \
	$(INSTALL_DATA) $(EXTRA_HEADERS) $(MTL_IMPS) $(includedir)


