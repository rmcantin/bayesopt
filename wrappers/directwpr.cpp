/*
-----------------------------------------------------------------------------
   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/
// BOOST Libraries
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include "directwpr.hpp"
#include "inneroptimization.hpp"
 
namespace DIRECT
{
  int evaluate_wrap_   (int *n, double *x, double *f, 
			int *flag__, int *iidata, 
			int *iisize, double *ddata, 
			int *idsize, char *cdata,
			int *icsize, int cdata_len)
  {
    array_adaptor<double> shared((*n), x);
    vector<double, array_adaptor<double> > sharedN((*n), shared); 
    
    // This is not very clever... but works!
    void *objPointer = iidata;
    InnerOptimization* OPTIMIZER = static_cast<InnerOptimization*>(objPointer);
    
    vector<double> vgrad(n);
    *f = OPTIMIZER->innerEvaluate(sharedN,vgrad);
    *flag__ = 0;
    
    return 0;
  } /* criteriawrap_ */

}
