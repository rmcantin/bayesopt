/** \mainpage Bayesian optimization toolbox
 * 
 * This is an efficient, C++ implementation of the Bayesian
 * optimization methodology for expensive functions, also called
 * Sequential Kriging Optimization (SKO), Efficient Global
 * Optimization (EGO). Basically, it uses a distribution over function
 * to build a metamodel of the unknown function for we are looking the
 * extrema, and then applying some sort of active learning strategy to
 * select the query points that provides most potential improvement in
 * the seek.
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
 * Originally, it was developed for the papers referenced at the end,
 * where a Gaussian process with hyperpriors on the mean and signal
 * covariance parameters. Then, the metamodel was constructed using
 * the Maximum a Posteriory (MAP) of the parameters.
 *
 * However, the library now has grown to support many more surrogate
 * models, with different distributions (Gaussian processes,
 * Student's-t processes, etc.), with many kernels and mean
 * functions. It also provides different criteria (even some combined
 * criteria) so the library can be used to any problem involving some
 * bounded optimization, stochastic bandits, active learning for
 * regression, etc.
 *
 * See more details at the
 * <a href="https://bitbucket.org/rmcantin/bayesian-optimization/wiki/Home">project
 * wiki</a>. 
 *
 * ----
 *
 * Ruben Martinez-Cantin, Nando de Freitas, Arnaud Doucet and Jose Castellanos.
 * Active Policy Learning for Robot Planning and Exploration under Uncertainty. 
 * Robotics: Science and Systems. 2007
 *
 * Ruben Martinez-Cantin, Nando de Freitas, Eric Brochu, Jose Castellanos and 
 * Arnaud Doucet (2009) A Bayesian Exploration-Exploitation Approach for Optimal
 * Online Sensing and Planning with a Visually Guided Mobile Robot. Autonomous 
 * Robots - Special Issue on Robot Learning, Part B, 27(3):93-103.
 *
 * ----
 * 
 * 
 * Copyright: See COPYING file that comes with this distribution

 */
