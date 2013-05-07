% Bayesian optimization using GPML and Finker's Direct
% 
% We assume that both files are loaded in the path
%
% Author: Ruben Martinez-Cantin <rmcantin@unizar.es>
%


clear all, close all
addpath('../testfunctions')
seed = 197; randn('seed',seed), rand('seed',seed)
plotting = 0;

f = @ackley;
bounds = ones(2,1)*[-32.768, 32.768];

% Set initial set of points using latin hypercube sampling
nInit = 30;   
nIter = 270;
xtr = lhsu(bounds(:,1),bounds(:,2),nInit);
ytr = zeros(size(xtr,1),1);
for i=1:size(xtr,1)
    ytr(i) = f(xtr(i,:));
end;

start_t = cputime;

% setup the GP
cov = {@covMaterniso,3}; 
sf = 1; ell = 0.4;  hyp0.cov  = log([ell;sf]);
mean = 'meanZero';       
sn = 0.005;          hyp0.lik  = log(sn);
lik = 'likGauss';     inf = 'infExact';

% Only learn the parameters once before doing optimization
Ncg = 50;
hyp = minimize(hyp0,'gp', -Ncg, inf, mean, cov, lik, xtr, ytr); % opt hypers


% set up DIRECT
opt.algorithm = NLOPT_GN_DIRECT_L;
opt.lower_bounds = bounds(:,1);
opt.upper_bounds = bounds(:,2);

opt.maxeval = 1000;
% opts.maxevals = 1000;
% opts.maxits = 1000;
% opts.showits = 0;
% Problem.f = criteria;

times = zeros(nIter,1);

for i = 1:nIter
    
    fprintf(1,'- Iter %i of %i \n', i,nIter);
    %[minval, newX, hist] = Direct(Problem, bounds, opts, hyp, inf, mean,...
    %                              cov, lik, xtr, ytr);
    criteria = (@(x) nei(x,hyp,inf,mean,cov,lik,xtr,ytr));
    opt.min_objective = criteria;
    
    [newX, minval, retcode] = nlopt_optimize(opt, [0 0]);
    newY = f(newX);
    xtr = [xtr; newX];    
    ytr = [ytr; newY];
    
    times(i) = cputime - start_t;
    
    fprintf(1,'New evaluation point. y = %2.4f \n', newY);
    disp(newX')

    if plotting
        % =====================================================================
        % PLOT RESULTS    =====================================================
        % =====================================================================
        x = 0:0.001:1;
        for j = 1:length(x)
            [crit(j), ymu(j), ys(j)] = feval(criteria, x(j) );
            yh(j) = f(x(j));
        end
    
        lw = 2;
        h = figure(1);
        clf;hold on;
    
        subplot(3,1,[1 2]);   hold on;
        title('Bayesian Optimization');
        plot(xtr,ytr,'ko','LineWidth',lw);
        plot(newX,newY,'ro','LineWidth',lw);
    
        plot(x,yh, 'b','LineWidth',lw);
        plot(x, ymu, 'k','LineWidth',lw);
        plot(x, ymu+2*ys, 'k:','LineWidth',lw);
        plot(x, ymu-2*ys, 'k:','LineWidth',lw);
        avalue = axis;   axis([0 1 avalue(3) avalue(4)]);
    
        subplot(3,1,3); hold on;
        title('Criteria');
        bar(x,-crit,'k');
        avalue = axis;   axis([0 1 avalue(3) avalue(4)]);
    
        pause(0.2);
    end
end

save('log_aiken.mat',times)

minY = min(ytr);





