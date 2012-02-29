 #ifndef INVERT_MATRIX_HPP
 #define INVERT_MATRIX_HPP

 // REMEMBER to update "lu.hpp" header includes from boost-CVS
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/vector_proxy.hpp>
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/triangular.hpp>
 #include <boost/numeric/ublas/lu.hpp>
 #include <boost/numeric/ublas/io.hpp>

 using namespace boost::numeric::ublas;
 
 /* Matrix inversion routine.
    Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class Min, class Mout>
int InvertMatrix (const Min& input, Mout& inverse) {
  
 	typedef permutation_matrix<std::size_t> pmatrix;
	typedef typename Min::value_type value_type;
 	// create a working copy of the input
 	matrix<value_type> A(input);
 	// create a permutation matrix for the LU-factorization
 	pmatrix pm(A.size1());

 	// perform LU-factorization
 	int res = lu_factorize(A,pm);
        if( res != 0 ) return res;

 	// create identity matrix of "inverse"
 	inverse.assign(identity_matrix<value_type>(A.size1()));

 	// backsubstitute to get the inverse
 	lu_substitute(A, pm, inverse);

 	return 0;
 }

 #endif //INVERT_MATRIX_HPP
