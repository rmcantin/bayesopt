#include <sstream>

#include "log.hpp"
#include "mean_functors.hpp"
#include "mean_atomic.hpp"
#include "mean_combined.hpp"



MeanFactory::MeanFactory()
{
  registry["mZero"] = & create_func<ZeroFunction>;
  registry["mOne"] = & create_func<OneFunction>;
  registry["mConst"] = & create_func<ConstantFunction>;
  registry["mLinear"] = & create_func<LinearFunction>;
  registry["mSum"] = & create_func<SumFunction>;
}


/** \brief Factory method for kernels. */
ParametricFunction* MeanFactory::create(mean_name name, 
					size_t input_dim)
{
  ParametricFunction* m_ptr;

  switch(name)
    {    
    case M_ZERO: m_ptr = new ZeroFunction(); break;
    case M_ONE: m_ptr = new OneFunction(); break;
    case M_CONSTANT: m_ptr = new ConstantFunction(); break;
    case M_LINEAR: m_ptr = new LinearFunction(); break;
    case M_LINEAR_CONSTANT: m_ptr = new LinearPlusConstantFunction(); break;
    default:
      FILE_LOG(logERROR) << "Error: mean function not supported.";
      return NULL;
    }

  m_ptr->init(input_dim);
  return m_ptr;
};

/** 
 * \brief Factory model for kernel functions
 * This function is based on the libgp library by Manuel Blum
 *      https://bitbucket.org/mblum/libgp
 * which follows the squeme of GPML by Rasmussen and Nickisch
 *     http://www.gaussianprocess.org/gpml/code/matlab/doc/
 * @param name string with the kernel structure
 * @param imput_dim number of input dimensions
 * @return kernel pointer
 */
ParametricFunction* MeanFactory::create(std::string name, size_t input_dim)
{
  ParametricFunction *mFunc;
  std::stringstream is(name);
  std::stringstream os(std::stringstream::out);
  std::stringstream os1(std::stringstream::out);
  std::stringstream os2(std::stringstream::out);
  char c;
  int i = 0, j = 0;
  while (is >> c) 
    {
      if (c == '(') i++;
      else if (c == ')') i--;
      else if (c == ',') j++;
      else 
	{
	  if (i == 0) os << c;
	  else if (j == 0) os1 << c;
	  else os2 << c;
	}
    }
  std::map<std::string,MeanFactory::create_func_definition>::iterator it = registry.find(os.str());
  if (it == registry.end()) 
    {
      FILE_LOG(logERROR) << "Error: Fatal error while parsing mean function: "
			 << os.str() << " not found" << std::endl;
      return NULL;
    } 
  mFunc = registry.find(os.str())->second();
  if (os1.str().length() == 0 && os2.str().length() == 0) 
    {
      mFunc->init(input_dim);
    } 
  else 
    {
      mFunc->init(input_dim, create(os1.str(),input_dim), 
		  create(os2.str(),input_dim));
    }
  return mFunc;

};



