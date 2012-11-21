/** \mainpage BayesOpt: A Bayesian optimization toolbox
 * 
 * This is an efficient, C++ implementation of the Bayesian
 * optimization methodology for nonlinear-optimization, experimental
 * design and stochastic bandits. In the literature it is also called
 * Sequential Kriging Optimization (SKO) or Efficient Global
 * Optimization (EGO). 
 * 
 * Basically, it uses a distribution over functions to build a
 * metamodel of the unknown function for we are looking the extrema,
 * and then apply some active learning strategy to select the query
 * points that provides most potential interest for the seek. For that
 * reason, it has been traditionally intended for optimization of
 * expensive function. However, the efficiency of the library make it
 * also interesting for many types of functions.
 *
 * It is intended to be both fast and clear for development and
 * research. At the same time, it does everything the "right way". For
 * example:
 * - latin hypercube sampling is use for the preliminary sampling step
 * - kernel parameters are trained with the preliminary samples and
 * fixed afterwards to avoid bias and divergence
 * - matrix algebra tricks are used to guarantee that any covariance
 * matrix remains SPD and reduce computational cost.
 * - etc.
 *
 * Originally, it was developed for as part of a robotics research
 * project \cite MartinezCantin09AR \cite MartinezCantin07RSS, where a
 * Gaussian process with hyperpriors on the mean and signal covariance
 * parameters. Then, the metamodel was constructed using the Maximum a
 * Posteriory (MAP) of the parameters.
 *
 * However, the library now has grown to support many more surrogate
 * models, with different distributions (Gaussian processes,
 * Student's-t processes, etc.), with many kernels and mean
 * functions. It also provides different criteria (even some combined
 * criteria) so the library can be used to any problem involving some
 * bounded optimization, stochastic bandits, active learning for
 * regression, etc.
 *
 * Start by reading the \ref install and the \ref reference
 *
 * You can also find more details at the
 * <a href="https://bitbucket.org/rmcantin/bayesian-optimization/wiki/Home">project
 * wiki</a>. 
 *
 * \b Important: This code is free to use. However, if you are using,
 * or plan to use, the library, specially if it is for research or
 * academic purposes, please send me an email at rmcantin@unizar.es
 * with your name, institution and a brief description of your
 * interest for this code (one or two lines).
 *
 * If you use BayesOpt in work that leads to a publication, we would
 * appreciate it if you would kindly cite BayesOpt in your
 * manuscript. Cite BayesOpt as something like:
 *
 * <hr>
 * Ruben Martinez-Cantin, <b>BayesOpt: a toolbox for nonlinear-optimization,
 * experimental design and stochastic bandits</b>,
 * http://bitbucket.org/rmcantin/bayesian-optimization
 * <hr>
 * 
 * Copyright: See COPYING file that comes with this distribution
 */
