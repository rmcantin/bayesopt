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

#ifndef  _KERNEL_FUNCTORS_HPP_
#define  _KERNEL_FUNCTORS_HPP_


#include "ctypes.h"
#include "specialtypes.hpp"
#include "elementwiseUblas.hpp"


////////////////////////////////////////////////////////////////////////////////
class Kernel
{
public:
  virtual void setScale( const vectord & theta ) {};
  virtual void setScale( double theta ) {};
  virtual double getScale(int index) {return 0.0;};
  virtual vectord getScale() {return zvectord(1);};
  virtual double operator()( const vectord &x1, const vectord &x2,
			     int grad_index = -1) = 0;
  virtual ~Kernel(){};
};


////////////////////////////////////////////////////////////////////////////////
class ISOkernel : public Kernel
{
public:
  void setScale( const vectord &theta ) {mTheta = theta(0);};
  void setScale( double theta ) {mTheta = theta;};
  double getScale(int index) {return mTheta;};
  vectord getScale() {return svectord(1,mTheta);};
  virtual ~ISOkernel(){};

protected:
  double mTheta;
};


////////////////////////////////////////////////////////////////////////////////
class ARDkernel : public Kernel
{
public:
  void setScale( const vectord &theta ){mTheta=theta;};
  double getScale(size_t index) {return mTheta(index);};
  vectord getScale() { return mTheta; };
  virtual ~ARDkernel(){};

protected:
  vectord mTheta;
};



////////////////////////////////////////////////////////////////////////////////
class MaternIso1: public ISOkernel
{
public:
  double operator()( const vectord &x1, const vectord &x2,
		     int grad_index = -1)
  {
    double r = norm_2(x1-x2)/mTheta;
    double er = exp(-r);

    if (grad_index < 0) 
      return er;
    else 
      return r*er;
  };
};


////////////////////////////////////////////////////////////////////////////////
class MaternIso3: public ISOkernel
{
public:
  double operator()( const vectord &x1, const vectord &x2,
		     int grad_index = -1 )
  {
    double r = sqrt(3) * norm_2(x1-x2)/mTheta;
    double er = exp(-r);

    if (grad_index < 0) 
      return (1+r)*er;
    else 
      return r*r*er; 
  };
};


////////////////////////////////////////////////////////////////////////////////
class MaternIso5: public ISOkernel
{
public:
  double operator()( const vectord &x1, const vectord &x2,
		     int grad_index = -1 )
  {
    double r = sqrt(5) * norm_2(x1-x2)/mTheta;
    double er = exp(-r);

    if (grad_index < 0) 
      return (1+r*(1+r/3))*er;
    else 
      return r*(1+r)/3*r*er; 
  };
};



////////////////////////////////////////////////////////////////////////////////
class SEIso: public ISOkernel
{
public:
  double operator()( const vectord &x1, const vectord &x2,
		     int grad_index = -1 )
  {
    double rl = norm_2(x1-x2)/mTheta;
    double k = rl*rl;

    if (grad_index < 0)
      return exp(-k/2);
    else 
      return exp(-k/2)*k;
  };
};


////////////////////////////////////////////////////////////////////////////////
class SEArd: public ARDkernel
{
public:
  double operator()( const vectord &x1, const vectord &x2,
		     int grad_index = -1 )
  {
    vectord xd = x1-x2;
    vectord ri = ublas_elementwise_div(xd, mTheta);

    double rl = norm_2(ri);
    double k = rl*rl;

    if (grad_index < 0)
      return exp(-k/2);
    else 
      return exp(-k/2)*sqrt(ri(grad_index));
  };
};

#endif
