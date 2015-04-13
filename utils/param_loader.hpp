
/**  \file param_loader.hpp \brief Allows to load parameters from file */
/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2014 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/


#ifndef  _PARAM_LOADER_HPP_
#define  _PARAM_LOADER_HPP_

#include "parameters.h"
#include "specialtypes.hpp"
#include "fileparser.hpp"


/**
 * Namespace of the library interface
 */
namespace bayesopt {

    //TODO (Javier): Complete class description
  /**
   * \brief .
   * 
   */
    class ParamLoader{
    public: 
        void loadFromFile(std::string filename, bopt_params &par);
        
        /* This one is useful to create param files from existing hard-coded params */
        void saveToFile(std::string filename, bopt_params &par);
    private:
        void loadOrSave(utils::FileParser &fp, bopt_params &par);
    };
} //namespace bayesopt


#endif
