


void MCMC::sliceSample(vectord &x)
{
  randFloat sample( mtRandom, realUniformDist(0,1) );
  size_t n = x.size();

  std::vector<int> perms = utils::returnIndexVector(n);
  utils::randomPerms(perms, mtRandom);

  for (size_t i = 0; i<n; ++i)
    {
      size_t ind = perms[i];
      double sigma = mSigma(ind);

      double y_max = obj->evaluate(x);
      double y = sample()*y_max;

      // Step out
      double x_cur = x(ind);
      double r = sample();
      double xl = x_cur - r * sigma;
      double xl = x_cur + (1-r)*sigma;

      if (mStepOut)
	{
	  x(ind) = xl;
	  while (obj->evaluate(x) > y) x(ind) -= sigma;
	  xl = x(ind);

	  x(ind) = xr;
	  while (obj->evaluate(x) > y) x(ind) += sigma;
	  xr = x(ind);
	}

      //Shrink
      bool on_slice = false;
      while (!on_slice)
	{
	  x(ind) = (xr-xl) * sample() + xl;
	  if (obj->evaluate(x) < y)
	    {
	      if      (x(ind) > x_cur)  xr = x(ind);
	      else if (x(ind) < x_cur)  xl = x(ind);
	      else throw std::runtime_error("Error in MCMC. Slice colapsed.");
	    }
	  else
	    {
	      on_slice = true;
	    }
	}
    }
}



