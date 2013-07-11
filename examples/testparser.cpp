#include <string>
#include <iostream>
#include "parser.hpp"

int main()
{
  std::string test = "One(Two, Three, Four(Five, Six))";
  std::string one;
  std::vector<std::string> vs;
  bayesopt::utils::parseExpresion(test,one,vs);

  std::cout << one << std::endl;
  for(std::vector<std::string>::iterator it = vs.begin();
      it != vs.end(); ++it)
    std::cout << *it << std::endl;
 
  return 0;
}
