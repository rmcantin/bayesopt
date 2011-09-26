/**/

int cubic_interpolation(double fa, double fb, 
			double xa, double xb, 
			double da, double db,
			double& root, double& A, double& B)
{
  A = 6*(fa-fb)+3*(db+da)*(xb-xa);
  B = 3*(fb-fa)-(2*da+db)*(xb-xa);
  root = B*B-A*da*(xb-xa);              // num. error possible, ok!
  if (root < 0.0)
    return -1;
  else
    {
      x3 = x1-d1*(x2-x1)^2/(B+sqrt(root));  
      return 1;
    }
}

#define INT  0.1 // don't reevaluate within 0.1 of the limit of current bracket
#define EXT  3.0           // extrapolate maximum 3 times the current step-size
#define MAX = 20                 // max 20 function evaluations per line search
#define RATIO = 10                               // maximum allowed slope ratio
#define SIG = 0.1
#define RHO = SIG/2    // SIG and RHO are the constants controlling the Wolfe-
// Powell conditions. SIG is the maximum allowed absolute ratio between
// previous and new slopes (derivatives in the search direction), thus setting
// SIG to low (positive) values forces higher precision in the line-searches.
// RHO is the minimum allowed fraction of the expected (from the slope at the
// initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
// Tuning of SIG (depending on the nature of the function to be optimized) may
// speed up the minimization; it is probably not worth playing much with RHO.

// The code falls naturally into 3 parts, after the initial line search is
// started in the direction of steepest descent. 1) we first enter a while loop
// which uses point 1 (p1) and (p2) to compute an extrapolation (p3), until we
// have extrapolated far enough (Wolfe-Powell conditions). 2) if necessary, we
// enter the second loop which takes p2, p3 and p4 chooses the subinterval
// containing a (local) minimum, and interpolates it, unil an acceptable point
// is found (Wolfe-Powell conditions). Note, that points are always maintained
// in order p0 <= p1 <= p2 < p3 < p4. 3) compute a new search direction using
// conjugate gradients (Polack-Ribiere flavour), or revert to steepest if there
// was a problem in the previous line-search. Return the best value so far, if
// two consecutive line-searches fail, or whenever we run out of function
// evaluations or line-searches. During extrapolation, the "f" function may fail
// either with an error or returning Nan or Inf, and minimize should handle this
// gracefully.

int minimize(vector<double>& X, int max_lines, int max_evals, 
	     int reduction = 1)
{
  int epoch = 0;                                 // zero the run length counters
  int evals = 0; 
  int ls_failed = 0;                       // no previous line search has failed

  double f0, d0;
  double x1, f1, d1;
  double x2, f2, d2;
  double x3, f3, d3;
  double x4, f4, d4;

  double A,B,root;

  vector<double> df0;
  maximum_likelihood(X,f0,df0);                // get function value and gradient

  epoch++;                                                      // count epochs?!
  vector<double> s = -df0;                 // initial search direction (steepest)
  d0 = dot(-s,s);                                                    // and slope
  x3 = reduction / (1-d0);                         // initial step is red/(|s|+1)

  while ((epoch < max_lines) && (evals < max_evals))
    {
      evals++;
      vector<double> X0 = X;                     // make a copy of current values
      double F0 = f0;
      vector<double> dF0 = df0;
      
      //TODO: Definir M

      while (1)                        // keep extrapolating as long as necessary
	{
	  x2 = 0; f2 = f0; d2 = d0; f3 = f0;
	  vector<double> df3 = df0;
	  
	  bool success = false;
	  while ((!success) && (M > 0))
	    {
	      M--; epoch++;                                     // count epochs?!
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
	  d3 = dot(df3,s);

	  if ((d3 > SIG*d0) || (f3 > f0+x3*RHO*d0) || (M == 0)) 
	    break;

	  x1 = x2; f1 = f2; d1 = d2;
	  x2 = x3; f2 = f3; d2 = d3;
	  A = 6*(f1-f2)+3*(d2+d1)*(x2-x1);
	  B = 3*(f2-f1)-(2*d1+d2)*(x2-x1);
	  root = B*B-A*d1*(x2-x1);                    // num. error possible, ok!
	  if (root < 0.0)
	    x3 = x2*EXT;                            // extrapolate maximum amount
	  else
	    {
	      x3 = x1-d1*(x2-x1)^2/(B+sqrt(root));    // num. error possible, ok!
	      if (x3 < 0) || (x3 > x2*EXT) // wrong sign or beyond extrap. limit?
		x3 = x2*EXT;                        // extrapolate maximum amount
	      else if x3 < x2+INT*(x2-x1)   // new point too close to prev point?
		x3 = x2+INT*(x2-x1);
	    }
	} /* while(1) */                                     // end extrapolation

      while (((abs(d3) > -SIG*d0) || (f3 > f0+x3*RHO*d0)) && (M > 0))
	{
	  if ( (d3 > 0) || (f3 > f0 + x3*RHO*d0) )
	    {
	      x4 = x3; f4 = f3; d4 = d3;
	    }
	  else 
	    {
	      x2 = x3; f2 = f3; d2 = d3;
	    }
	  if (f4 > f0)   // TODO: Check bug. f4 might be undefined.
	    {
	      A = (f4-f2-d2*(x4-x2));         // quad interpolation
	      if (A == 0.0)
		x3 = (x2+x4)/2;     // if we had a numerical problem then bisect
	      else
		x3 = x2-(0.5*d2*(x4-x2)^2)/A;
	    }
	  else
	    {
	      A = 6*(f2-f4)/(x4-x2)+3*(d4+d2);            // cubic interpolation
	      B = 3*(f4-f2)-(2*d2+d4)*(x4-x2);
	      root = B*B-A*d1*(x2-x1);               // num. error possible, ok!
	      if ((root < 0.0) || (A == 0.0))
		x3 = (x2+x4)/2;     // if we had a numerical problem then bisect
	      else
		x3 = x2+(sqrt(B*B-A*d2*(x4-x2)^2)-B)/A;     
	    }
	  x3 = max(min(x3,x4-INT*(x4-x2)),x2+INT*(s4-x2));
	  maximum_likelihood(X+x3*s,f3,df3); // get function value and gradient
	  if (f3 < F0)
	    {
	      X0 = X+x3*s;
	      F0 = f3;
	      dF0 = df3;
	    }
	  M--; epochs++;
	  d3 = dot(df3,s);
	}

      if ((fabs(d3) < -SIG*d0) && (f3 < f0+x3+RHO+d0)) //if line search succeeded
	{
	  X = X+dot(x3,s); f0 = f3;
	  s = dot(df3,df3) - dot(df0,df3) 
	}

    }  /* while(epoch... */


