/**/

#define INT  0.1    // don't reevaluate within 0.1 of the limit of the current bracket
#define EXT  3.0                  // extrapolate maximum 3 times the current step-size
#define MAX = 20                         // max 20 function evaluations per line search
#define RATIO = 10                                       // maximum allowed slope ratio
#define SIG = 0.1
#define RHO = SIG/2    // SIG and RHO are the constants controlling the Wolfe-

int minimize(vector<double>& result, fpointer f, int length, int red = 1)
{
  int i = 0;
  int ls_failed = 0;
  
  vector<double> f0;
  matrix<double> df0;

  f(result,f0,df0)


