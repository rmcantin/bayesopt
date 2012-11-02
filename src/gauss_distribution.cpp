#include <boost/math/special_functions/factorials.hpp>
#include "gauss_distribution.hpp"

GaussianDistribution::GaussianDistribution()
{
  mean_ = 0.0;  std_ = 1.0;
}


GaussianDistribution::~GaussianDistribution(){}

double GaussianDistribution::negativeExpectedImprovement(double min,
							 size_t g)
{
  
  using boost::math::factorial;

  const double diff = min - mean_;
  const double z = diff / std_;
  const double pdf_z = pdf(d_,z);
  const double cdf_z = cdf(d_,z);
  
  if (g == 1)
    return -1.0 * ( diff * cdf_z + std_ * pdf_z );
  else
    {
      const double fg = factorial<double>(g);

      double Tm2 = cdf_z;
      double Tm1 = pdf_z;
     
      double Tact;
      double sumEI = pow(z,g)*Tm2 - g*pow(z,g-1)*Tm1;

      for (size_t ii = 2; ii < g; ++ii) 
	{
	  Tact = (ii-1)*Tm2 - pdf_z*pow(z,ii-1);
	  sumEI += pow(-1.0,ii)* 
	    (fg / ( factorial<double>(ii)*factorial<double>(g-ii) ) )*
	    pow(z,g-ii)*Tact;
	  
	  //roll-up
	  Tm2 = Tm1;   Tm1 = Tact;
	}
      return -1.0 * pow(std_,g) * sumEI;
    }
  
}  // negativeExpectedImprovement

double GaussianDistribution::lowerConfidenceBound(double beta)
{    
  return mean_ - beta*std_;
}  // lowerConfidenceBound


double GaussianDistribution::negativeProbabilityOfImprovement(double min,
							      double epsilon)
{
  return -cdf(d_,(min - mean_ + epsilon)/std_);
}  // negativeProbabilityOfImprovement


double GaussianDistribution::sample_query(randEngine& eng)
{ 
  randNFloat sample(eng,normalDist(mean_,std_));
  return sample();
} // sample_query
