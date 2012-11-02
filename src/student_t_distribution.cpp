#include "student_t_distribution.hpp"

StudentTDistribution::StudentTDistribution(size_t dof):
  d_(dof)
{
  mean_ = 0.0;  std_ = 1.0; dof_ = dof;
}

StudentTDistribution::~StudentTDistribution(){}


double StudentTDistribution::negativeExpectedImprovement(double min,
							 size_t g)
{
  const double diff = min - mean_;
  const double z = diff / std_;
  
  assert((g == 1) && "Students t EI with exponent not yet supported.");
  return -(diff * cdf(d_,z) +  (dof_*std_+z*diff)/(dof_-1) * pdf(d_,z) ); 
}  // negativeExpectedImprovement


double StudentTDistribution::lowerConfidenceBound(double beta)
{    
  return mean_ - beta*std_/sqrt(dof_);
}  // lowerConfidenceBound

double StudentTDistribution::negativeProbabilityOfImprovement(double min,
							      double epsilon)
{  
  return -cdf(d_,(min - mean_ + epsilon)/std_);
}  // negativeProbabilityOfImprovement


double StudentTDistribution::sample_query(randEngine& eng)
{ 
  double n = static_cast<double>(dof_);
  randNFloat normal(eng,normalDist(mean_,std_));
  randGFloat gamma(eng,gammaDist(n/2.0));
  return normal() / sqrt(2*gamma()/n);
}  // sample_query
