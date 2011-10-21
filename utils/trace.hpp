#ifndef __TRACE_HPP__
#define __TRACE_HPP__

//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "specialtypes.hpp"

template<class E>
typename E::value_type trace(const E &A)
{
  size_t n = std::min(A.size1(),A.size2());
  typename E::value_type sum = 0;
  for (size_t i=0; i<n; ++i)
    sum += A(i,i);

  return sum; 
}

template<class E1, class E2>
typename E1::value_type trace_prod(const E1 & A, const E2 & B )
{
  size_t n = std::min(A.size1(),B.size2());
  typename E1::value_type sum = 0;
  for (size_t i=0; i<n; ++i)
    sum += ublas::inner_prod(ublas::row(A,i),ublas::column(B,i));

  return sum; 
}



#endif
