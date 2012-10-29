#include "mean_functors.hpp"

ParametricFunction* ParametricFunction::create(mean_name name, const vectord& params)
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
      std::cout << "Error: mean function not supported." << std::endl;
      return NULL;
    }

  m_ptr->setParameters(params);
  return m_ptr;
};

