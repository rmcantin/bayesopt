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

#include "bopt_state.hpp"
#include "fileparser.hpp"

#include <fstream>
#include <iostream>

namespace bayesopt
{       
    void BOptState::saveToFile(std::string filename){
        loadOrSave(filename, false);
    }
    
    void BOptState::loadFromFile(std::string filename){
        loadOrSave(filename, true);
        
        std::cout << "mCurrentIter=" << mCurrentIter << std::endl;
        std::cout << "mCounterStuck=" << mCounterStuck << std::endl;
        std::cout << "mYPrev=" << std::setprecision(10) << mYPrev << std::endl;
        std::cout << "mParameters.kernel.n_hp=" << mParameters.kernel.n_hp << std::endl;
    }
    
    void BOptState::loadOrSave(std::string filename, bool readMode){
        utils::FileParser fp(filename, readMode);        
        fp.readOrWrite("mCurrentIter", mCurrentIter, readMode);
        fp.readOrWrite("mCounterStuck", mCounterStuck, readMode);
        fp.readOrWrite("mYPrev", mYPrev, readMode);
        fp.readOrWrite("mParameters.kernel.n_hp", mParameters.kernel.n_hp, readMode);
    } 
} //namespace bayesopt

