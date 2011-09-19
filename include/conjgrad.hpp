/**/

#define INT  0.1    // don't reevaluate within 0.1 of the limit of the current bracket
#define EXT  3.0                  // extrapolate maximum 3 times the current step-size
#define MAX = 20                         // max 20 function evaluations per line search
#define RATIO = 10                                       // maximum allowed slope ratio
#define SIG = 0.1
#define RHO = SIG/2    // SIG and RHO are the constants controlling the Wolfe-

int minimize(vector<double>& result, int max_lines, int max_evals, 
	     int reduction = 1)
{
  int epoch = 0;
  int evals = 0;
  int ls_failed = 0;

  double f0;

  vector<double> df0;
  
  maximum_likelihood(result,f0,df0);

  epoch++;
  vector<double> s = -df0;
  double d0 = dot(-s,s);
  double x3 = reduction / (1-d0);

  while ((epoch < max_lines) && (evals < max_evals))
    {
      evals++;
      vector<double> X0 = result;
      double F0 = f0;
      vector<double> dF0 = df0;
      
      //TODO: Definir M

      while (1)
	{
	  double x2 = 0;
	  double f2 = f0;
	  double d2 = d0;
	  double f3 = f0;
	  vector<double> df3 = df0;
	  
	  bool success = false;
	  while ((!success) && (M > 0))
	    {
	      M--; epoch++;
	      if (maximum_likelihood(result+x3,f3,df3) == 1)
		success = true;
	      else x3 = (x2+x3)/2;
	    }
	  if (f3 < F0)
	    {
	      X0 = X+x3*s;
	      F0 = f3;
	      dF0 = df3;
	    }
	}
    }


