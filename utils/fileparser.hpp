/**  \file fileparser.hpp \brief Functions to write and parse data files */
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

#ifndef  _FILEPARSER_HPP_
#define  _FILEPARSER_HPP_

#include <vector>

#include <fstream>
#include <iostream>

#include <sstream>

#include <iomanip>

namespace bayesopt 
{  
  namespace utils 
  {
    class FileParser{
    public:
        FileParser(std::string filename);
        FileParser(std::string filename, bool readMode);
        ~FileParser();
        
        /* Write stream and read stream open/close functions*/
        void openOutput();
        void openInput();
        void close();
        
        /* Data write/read function*/
        void write(std::string name, std::string value);
        void read(std::string name, std::string &value);
        std::string read(std::string name);
        
        /* Array write/read function */
        void write(std::string name, const std::vector<std::string> &arr, const std::vector<int> &dims);
        void read(std::string name, std::vector<std::string> &arr, std::vector<int> &dims);
        
        /* Template definitions and implementation */
        template <typename T> void write(std::string name, T value){
            std::string str = FileParser::to_string<T>(value);
            write(name, str);
        }
        template <typename T> void read(std::string name, T &value){
            value = to_number<T>(FileParser::read(name));
        }
        template <typename T> void readOrWrite(std::string name, T&value, bool readMode){
            if(readMode){
                read(name, value);
            }
            else{
                T v = value;
                write(name,v);
            }
        }
        
    private:
        /* Search variable in file */
        bool movePointer(std::string name, std::string &content);
        bool startsWith(std::string all, std::string sub);
        
        /* Parses the contents as an array */
        void parseArray(std::string contents, std::vector<std::string> &arr, std::vector<int> &dims);
        
        /* private members */
        std::string filename;
        std::ofstream output;
        std::ifstream input;
        
        std::string currentLine;
        
        /* Template definitions and implementation */
        // TODO (Javier): std to_string, stoi, etc... not working on MinGW and does not provide precision, workaround:
        template <typename T>
        std::string to_string(T value, int prec = 10)
        {
          std::ostringstream os ;
          os << std::setprecision(prec) << value ;
          return os.str() ;
        }
        
        template <typename T>
        T to_number(std::string str, int prec = 10)
        {
            std::istringstream ss(str);
            T result;
            return ss >> std::setprecision(prec) >> result ? result : 0;
        }
    };
  } //namespace utils
} //namespace bayesopt

#endif
