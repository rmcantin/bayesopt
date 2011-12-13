#include <valarray>
#include "bayesopt.hpp"

inline double sqr( double x ){ return x*x; }

class TestOneD: public SKO
{
public:
  TestOneD(NonParametricProcess* gp):
  SKO(gp) {}

  double evaluateSample( const vectord& xin)
  {
    if (xin.size() > 1)
      {
	std::cout << "This only works for 1D inputs" << std::endl;
	return 0.0;
      }

    double x = xin(0);
    return sqr(x-0.5) + sin(20*x)*0.2;
  };

  bool checkReachability(const vectord &query)
  {return true;};

};

int main(int nargs, char *args[])
{
  size_t nIterations = 300;
  NonParametricProcess* gp = new BasicGaussianProcess(1.0,DEF_REGULARIZER);
  TestOneD opt(gp);
  vectord result(1);
  opt.setCriteria(c_ei);

  opt.optimize(result,nIterations);
  
  return 1;
}
