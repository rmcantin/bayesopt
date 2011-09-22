//
// C++ Implementation: meanfuncs
//
// Description: 
//
//
// Author: Ruben Martinez-Cantin  <rmcantin@unizar.es>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//



#ifndef  _MEANFUNCS_HPP_
#define  _MEANFUNCS_HPP_

#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;	

namespace means
{
  double One( const vector<double> &x )
  { return 1; } 
  
  double Zero( const vector<double> &x)
  { return 0; } 

  double Linear (const vector<double> &x,
		 const vector<double> &a)
  { return inner_prod(x,a); }
}


#endif
