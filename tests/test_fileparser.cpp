#include <iostream>
#include "fileparser.hpp"
#include "ublas_extra.hpp"

int main()
{
    std::string filename = "test_fileparser.dat";
    
    size_t var_size_t = 24;
    int var_int_signed = -32;   // This 2 varaibles are necessary to check that FileParser
    int var_int = 25;           // does not retrieve "int_signed" when searching for only "int"
    double var_double = 1.2345678;
    char var_char_arr[5] = {'c', 'h', 'a', 'r'};
    std::string var_string = "string_value";
    double var_double_arr[5] = {0,1.1,3.333,5.55555,8.88888888};
    
    
    bayesopt::utils::FileParser fp(filename);
    // Write values
    fp.openOutput();
    
    fp.write("size_t", var_size_t);
    fp.write("int_signed", var_int_signed);
    fp.write("int", var_int);
    fp.write("double", var_double);
    fp.write_chars("char_arr", var_char_arr);
    fp.write("string", var_string);
    fp.write_double_array("double_arr", var_double_arr,5);    
    
    fp.close();
    
    // Create new variables
    size_t cp_size_t;
    int cp_int_signed;
    int cp_int;
    double cp_double;
    char cp_char_arr[5];
    std::string cp_string;
    double cp_double_arr[5];
    
    // Read values
    fp.openInput();
    
    fp.read("size_t", cp_size_t);
    fp.read("int_signed", cp_int_signed);
    fp.read("int", cp_int);
    fp.read("double", cp_double);
    fp.read_chars("char_arr", cp_char_arr);
    fp.read("string", cp_string);
    fp.read_double_array("double_arr", cp_double_arr,5);
    
    fp.close();
    
    // Display .dat file
    std::cout << "Displaying saved file contents:" << std::endl;
    std::ifstream input;
    input.open(filename.c_str());
    std::string line;
    while(getline(input,line)){
        std::cout << line << std::endl;
    }
    input.close();
    std::cout << std::endl;
    
    // Display variable contents
    std::cout.precision(10);
    std::cout << "Displaying variable contents:" << std::endl;
    std::cout << "size_t=" << cp_size_t << std::endl;
    std::cout << "int_signed=" << cp_int_signed << std::endl;
    std::cout << "int=" << cp_int << std::endl;
    std::cout << "double=" << cp_double << std::endl;
    std::cout << "char_arr=" << std::string(cp_char_arr) << std::endl;
    std::cout << "string=" << cp_string << std::endl;
    std::cout << "double_arr=" << bayesopt::utils::array2vector(cp_double_arr, 5) << std::endl;
    std::cout << std::endl;
    
    // Try to remove used .dat file
    if( remove( filename.c_str() ) == 0 ){
        std::cout << "File \"" << filename << "\" successfully file" << std::endl;
    }
    else{
        std::cout << "Error: cannot remove \"" << filename << "\" file" << std::endl; 
    }
    
    // Some Values Checks
    int returnValue = 0;
    if(var_size_t != cp_size_t){
        returnValue = -1;
        std::cout << "ERROR: expected " << var_size_t << " to be equals to " << cp_size_t << std::endl;
    }
    if(var_int_signed != cp_int_signed){
        returnValue = -1;
        std::cout << "ERROR: expected " << var_int_signed << " to be equals to " << cp_int_signed << std::endl;
    }
    if(var_int != cp_int){
        returnValue = -1;
        std::cout << "ERROR: expected " << var_int << " to be equals to " << cp_int << std::endl;
    }
    if(var_double != cp_double){
        returnValue = -1;
        std::cout << "ERROR: expected " << var_double << " to be equals to " << cp_double << std::endl;
        std::cout << "(Note: check values, could be a precision error)" << std::endl;
    }
    if(var_string != cp_string){
        returnValue = -1;
        std::cout << "ERROR: expected " << var_string << " to be equals to " << cp_string << std::endl;
    }
    return returnValue;
}
