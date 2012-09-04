/*
-----------------------------------------------------------------------------
   This file is part of BayesOptimization, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOptimization is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOptimization is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOptimization.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/

#ifndef __LOGGER_HPP__
#define __LOGGER_HPP__

#include <iostream>
#include <fstream>

class Logger
{
public:
  Logger( bool useFile = false, 
	 const char* filename = "bopt_output.log"):
    mOutput(useFile ? mFileOutput : std::cout)
  {
    if (useFile)
    {
      mFileOutput.open(filename);
      if ( !mFileOutput.is_open() ) 
	{
	  std::cout << "Cannot open log file." << std::endl;
	}
    }

  }

  virtual ~Logger()
  {
    if ( mFileOutput.is_open() )
      mFileOutput.close();
  };

protected:
  std::ofstream mFileOutput;
  std::ostream& mOutput;           ///< Log file in case it is used
};

#endif
