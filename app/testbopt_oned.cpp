#include <valarray>
#include "bayesoptcont.hpp"

inline double sqr( double x ){ return x*x; }

class TestOneD: public BayesOptContinuous
{
public:
  TestOneD(bopt_params par):
  BayesOptContinuous(par) {}

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
  bopt_params parameters = initialize_parameters_to_default();
  parameters.n_iterations = 300;
  parameters.theta[0] = 1.0;
  parameters.n_theta = 1;
  parameters.c_name = C_GP_HEDGE;
  TestOneD opt(parameters);
  vectord result(1);
  opt.optimize(result);
  
  return 1;
}
