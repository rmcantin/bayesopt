
#ifndef  _RANDGEN_HPP_
#define  _RANDGEN_HPP_

#include <boost/version.hpp>
#include <boost/random.hpp>

// Types for pseudorandom number generators.
// According to documentation, the uniform floating point distribution is buggy.
// I've read lots of webpages complaining, but everyone uses it.

typedef boost::mt19937                                              randEngine;

typedef boost::uniform_real<>				       realUniformDist;
typedef boost::uniform_int<> 					intUniformDist;
typedef boost::normal_distribution<>                                normalDist;
typedef boost::gamma_distribution<>                                  gammaDist;

// #if BOOST_VERSION >= 104700
//    typedef boost::student_t_distribution<>                           studTDist;

typedef boost::variate_generator<randEngine&, intUniformDist>          randInt;
typedef boost::variate_generator<randEngine&, normalDist>           randNFloat;
typedef boost::variate_generator<randEngine&, gammaDist>            randGFloat;
typedef boost::variate_generator<randEngine&, realUniformDist>       randFloat;

#endif
