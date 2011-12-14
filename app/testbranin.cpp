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
    if (xin.size() != 2)
      {
	std::cout << "This only works for 2D inputs" << std::endl;
	return 0.0;
      }

    double x = xin(0) * 15 - 5;
    double y = xin(1) * 15;

    return sqr(y-(5.1/(4*sqr(M_PI)))*sqr(x)+5*x/M_PI-6)+10*(1-1/(8*M_PI))*cos(x)+10;
  };

  bool checkReachability(const vectord &query)
  {return true;};

};

int main(int nargs, char *args[])
{
  size_t nIterations = 400;
  NonParametricProcess* gp = new StudentTProcess(1.0,DEF_REGULARIZER);
  TestOneD opt(gp);
  vectord result(2);
  opt.setCriteria(c_gp_hedge);

  opt.optimize(result,nIterations);

  std::cout << "Soluciones: 0.1239    0.8183" << std::endl;
  std::cout << "Soluciones: 0.1239    0.1517" << std::endl;
  std::cout << "Soluciones: 0.9617    0.1650" << std::endl;
  
  return 1;
}
