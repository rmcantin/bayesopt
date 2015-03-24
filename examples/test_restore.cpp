#include "testfunctions.hpp"
#include "bopt_state.hpp"

int main()
{
    std::cout << "RESTORING OPTIMIZATION" << std::endl;
    std::cout << "======================" << std::endl;    
    
    // Second optimization (restored from first optimization state)
    bopt_params par2 = initialize_parameters_to_default();
    par2.n_iterations = 190;
    par2.random_seed = 0;
    par2.verbose_level = 1;
    par2.noise = 1e-10;

    BraninNormalized branin2(par2);

    // Restore operation and optimization
    bayesopt::BOptState state2;
    state2.loadFromFile("state.dat");
    branin2.restoreOptimization(state2);
    for(size_t i = branin2.getCurrentIter(); i < par2.n_iterations; i++){
        branin2.stepOptimization();
    }
    
    vectord result = branin2.getFinalResult();
    std::cout << "Branin2 Result: " << result << "->" 
        << branin2.evaluateSample(result) << std::endl;
}
