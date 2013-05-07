function [ei, varargout] = nei(query, hyp, inf, mean, cov,...
                               lik, xtr, ytr, params )
% NEGEXPIMPROVEMENT Compute the *negative* expected improvement.
%
% DIRECT will call this with the query vector as a column, but we want
% it to be a row.  Make sure we can call it both ways.

if (size(query,1) > size(query,2))
   query = query.';
end

[ymu, ys2] = gp(hyp, inf, mean, cov, lik, xtr, ytr, query);  % predict


s = sqrt(ys2);             % Predicted std
yres = min(ytr) - ymu;     % Predicted error
ynorm = yres/s;            % Normalize error

ei = yres * normcdf(ynorm) + s * normpdf(ynorm);
ei = max([0 ei]);

% DIRECT search for a minimum, and we want the maximum expected
% improvement, so...
ei = -ei;

varargout{1} = ymu;    varargout{2} = s;