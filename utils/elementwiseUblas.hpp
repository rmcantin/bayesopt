#ifndef  _ELEMENTWISE_UBLAS_HPP_
#define  _ELEMENTWISE_UBLAS_HPP_

// BOOST Libraries
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <algorithm>


template<class type_data>
type_data op_sum (type_data i, type_data j) { return i+j; }

template<class type_data>
type_data op_prod (type_data i, type_data j) { return i*j; }

template <class T>
T ublas_elementwise_add(T a, T b)
{
  T c(a);

  for(size_t ii = 0; ii < c.size(); ii++)
      c(ii) += b(ii);
    
  return c;
 
}

template <class T>
T ublas_elementwise_prod(T a, T b)
{
  T c(a);

  for(size_t ii = 0; ii < c.size(); ii++)
      c(ii) *= b(ii);
    
  return c;
 
}


#endif
