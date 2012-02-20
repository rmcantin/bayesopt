% Usage: [xmin, fmin] = bayesopt(@function_handler, nDimensions, params)
%        [xmin, fmin] = bayesopt(@function_handler, nDimensions, params,
%                                lowerBound, upperBound) 
%
%        [xmin, fmin] = bayesopt('function_name', nDimensions, params)
%        [xmin, fmin] = bayesopt('function_name', nDimensions, params,
%                                lowerBound, upperBound) 
%
%
% nDimensions is the number of dimensions of the query vector.
%
% Params is a struct which may have the fields 
%    theta : kernel parameters
%    alpha, beta : std hyperprior
%    delta : mean hyperprior
%    noise : likelihood std
%    iterations : maximum number of iterations
%
% lowerBound and upperBound should be a nDim x 1 or 1 x nDim vectors with
%      the lower and upper bound for each component. (optional, default 0-1)
%
