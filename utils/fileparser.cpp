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

#include "fileparser.hpp"

namespace bayesopt
{   
  namespace utils 
  {
    FileParser::FileParser(std::string filename)
    : filename(filename), input(), output(){}
    
    FileParser::FileParser(std::string filename, bool readMode)
    : filename(filename), input(), output(){
        if(readMode){
            openInput();
        }
        else{
            openOutput();
        }
    }

    FileParser::~FileParser(){
        close();
    }
    
    void FileParser::openOutput(){
        close();
        output.open(filename.c_str());
    }
    void FileParser::openInput(){
        close();
        input.open(filename.c_str());
    }
    void FileParser::close(){
        output.close();
        input.close();
        currentLine = "";
    }
    
    void FileParser::write(std::string name, std::string value){
        output << name << "=" << value << std::endl;
    }
    void FileParser::read(std::string name, std::string &value){
        if(!movePointer(name,value)){
            std::cerr << "Variable: " << name << " does not exist in file: " << filename << std::endl;
        }
    }
    std::string FileParser::read(std::string name){
        std::string ret;
        read(name, ret);
        return ret;
    }
        
    void FileParser::write(std::string name, const std::vector<std::string> &arr, const std::vector<int> &dims){
        // Write dimensions
        output << name << "=(";
        for(std::vector<int>::const_iterator it = dims.begin(); it != dims.end(); ++it) {
            if(it != dims.begin()){
                output << " ";
            }
            output << *it;
        }
        output << ")";
        
        // Write array
        for(std::vector<std::string>::const_iterator it = arr.begin(); it != arr.end(); ++it) {
            if(it != arr.begin()){
                output << " ";
            }
            output << *it;
        }
        output << std::endl;
    }
    void FileParser::read(std::string name, std::vector<std::string> &arr, std::vector<int> &dims){
        std::string contents;
        if(movePointer(name,contents)){
            parseArray(contents, arr, dims);
        }
        else{
            std::cerr << "Variable: " << name << " does not exist in file: " << filename << std::endl;
        }
    }
    
    bool FileParser::movePointer(std::string name, std::string &contents){
        std::cout << "DEBUG::movePointer| variable name: " << name << std::endl;
        
        std::cout << "DEBUG::movePointer| currentLine: " << currentLine << std::endl;
        if(currentLine.length() > 0 && startsWith(currentLine, name)){
            contents = currentLine.substr(name.length()+1);
            std::cout << "DEBUG::movePointer| returnvalue" << contents << std::endl;
            return true;
        }
        
        int debug_iters = 0;
        
        // Wrap the file around in the search of variable
        for(int i=0; i<2; i++){
            while (getline( input, currentLine)){
                if(currentLine.length() > 0 && startsWith(currentLine, name)){
                    contents = currentLine.substr(name.length()+1);
                    std::cout << "DEBUG::movePointer| returnvalue:" << contents << std::endl;
                    std::cout << "DEBUG::movePointer| totalIters:" << debug_iters << std::endl;
                    return true;
                }
                debug_iters++;
            }
            input.clear();
            input.seekg(0, std::ios::beg);
            std::cout << "DEBUG::searchVariable| EOF detected" << std::endl;
            std::cout << "DEBUG::movePointer| totalIters:" << debug_iters << std::endl;
        }
        contents = "";
        return false;
        
    }
    
    bool FileParser::startsWith(std::string all, std::string sub){
        return all.rfind(sub, 0) == 0;
    }
    
    void FileParser::parseArray(std::string contents, std::vector<std::string> &arr, std::vector<int> &dims){
        std::cout << "DEBUG::parseArray| not implemented" << std::endl;
    }
  } //namespace utils
} //namespace bayesopt

