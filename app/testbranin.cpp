#include <valarray>
#include "bayesoptcont.hpp"

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

inline double sqr( double x ){ return x*x; }

class TestOneD: public BayesOptContinuous
{
public:
  TestOneD(bopt_params par):
  BayesOptContinuous(par) {}

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
  bopt_params par = initialize_parameters_to_default();
  par.n_iterations = 200;
  par.n_init_samples = 50;
  par.theta[0] = 1.0;
  par.n_theta = 1;
  par.c_name = C_GP_HEDGE;
  TestOneD opt(par);
  vectord result(2);

  opt.optimize(result);

  std::cout << "Soluciones: 0.1239    0.8183" << std::endl;
  std::cout << "Soluciones: 0.1239    0.1517" << std::endl;
  std::cout << "Soluciones: 0.9617    0.1650" << std::endl;
  
  return 1;
}
