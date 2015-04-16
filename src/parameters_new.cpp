/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2014 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/

#include "parameters_new.hpp"

namespace bayesopt
{       
    /*-----------------------------------------------------------*/
    /* Default parameters                                        */
    /*-----------------------------------------------------------*/

    /* Nonparametric process "parameters" */
    const std::string KERNEL_NAME = "kMaternARD5";
    const double KERNEL_THETA    = 1.0;
    const double KERNEL_SIGMA    = 10.0;
    const std::string MEAN_NAME = "mConst";
    const double MEAN_MU         = 1.0;
    const double MEAN_SIGMA      = 1000.0;
    const double PRIOR_ALPHA     = 1.0;
    const double PRIOR_BETA      = 1.0;
    const double DEFAULT_SIGMA   = 1.0;
    const double DEFAULT_NOISE   = 1e-6;

    /* Algorithm parameters */
    const size_t DEFAULT_ITERATIONS         = 190;
    const size_t DEFAULT_INNER_EVALUATIONS  = 500; /**< Used per dimmension */
    const size_t DEFAULT_INIT_SAMPLES       = 10;
    const size_t DEFAULT_ITERATIONS_RELEARN = 50;
    
    /* Logging and Files */
    const size_t DEFAULT_VERBOSE    = 1;
    const std::string LOG_FILENAME  = "bayesopt.log";
    const std::string SAVE_FILENAME = "bayesopt.dat";
    const std::string LOAD_FILENAME = "bayesopt.dat";
    
    const std::string SURR_NAME = "sGaussianProcess";
    const std::string CRIT_NAME = "cEI";
    
    /*
     * KernelParameters Class
     */
    KernelParameters::KernelParameters():
        hp_mean(), hp_std(){
        // Set default values
        name = KERNEL_NAME;
        hp_mean.push_back(KERNEL_THETA);
        hp_std.push_back(KERNEL_SIGMA);
        n_hp = 1;
    }
    
    /*
     * MeanParameters Class
     */
    MeanParameters::MeanParameters():
    coef_mean(), coef_std(){
        // Set default values
        name = MEAN_NAME;
        coef_mean.push_back(MEAN_MU);
        coef_std.push_back(MEAN_SIGMA);
        n_coef = 1;
    }
    
    /*
     * Parameters Class
     */
    Parameters::Parameters():
        kernel(), mean(){
        // Set default values
        init_default();
    }
        
    Parameters::Parameters(bopt_params c_params):
        kernel(), mean(){
        //TODO (Javier): get values from the bopt_params struct
    }
    
    bopt_params Parameters::generate_bopt_params(){
        bopt_params c_params;
        //TODO (Javier): fill bopt_params struct values
        return c_params;
    }
    
    void Parameters::set_learning(std::string name){
        l_type = str2learn(name.c_str());
    }
    std::string Parameters::get_learning(){
        return std::string(learn2str(l_type));
    }

    void Parameters::set_score(std::string name){
        sc_type = str2score(name.c_str());
    }
    std::string Parameters::get_score(){
        return std::string(score2str(sc_type));
    }

    void Parameters::init_default(){
        n_iterations = DEFAULT_ITERATIONS;
        n_inner_iterations = DEFAULT_INNER_EVALUATIONS;
        n_init_samples = DEFAULT_INIT_SAMPLES;
        n_iter_relearn = DEFAULT_ITERATIONS_RELEARN;
        
        init_method = 1;
        random_seed = -1;
        
        verbose_level = DEFAULT_VERBOSE;
        log_filename = LOG_FILENAME;
        
        load_save_flag = 0;
        load_filename = LOAD_FILENAME;
        save_filename = SAVE_FILENAME;
        
        surr_name = SURR_NAME;
        
        sigma_s = DEFAULT_SIGMA;
        noise = DEFAULT_NOISE;
        alpha = PRIOR_ALPHA;
        beta = PRIOR_BETA;
        
        l_all = false;
        l_type = L_EMPIRICAL;
        sc_type = SC_MAP;
        
        epsilon = 0.0;
        force_jump = 20;
        
        crit_name = CRIT_NAME;
        n_crit_params = 0;
    }
} //namespace bayesopt

