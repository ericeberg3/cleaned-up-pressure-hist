function  [x_keep, L_keep, count] = mcmc(func,data,x0,xstep,xbnds,sigma,Niter,varargin)
%
% [x_keep, L_keep, count] = mcmc(func,data,x0,xstep,sigma,Niter,varargin)
%
% subroutine for MCMC sampling using Metropolis-Hasting w/ normal
% distribution.
%
% Inputs:
%       func = function name eg:   mcmc('travel_time',....)
%               parameter required for function in varargin
%       data = vector of observations
%       x0   = initial estimate of parameter vector
%       xstep = step size in all parameter directions
%       xbnds = bounds (
%       sigma = sigma of normal distribution 
%       Niter = number of iterations
%
% Outputs:
%       x_keep = array of samples
%       L_keep = likelihood of samples
%       count  = number of accepted. Acceptance ratio is count/Niter
%
%  P Segall: 2012
% TODO? Should I add something to seed rand? Better if this is done in
% calling program
% rng('shuffle');

fun = fcnchk(func);
% fun = str2func(func);

%number of elements in x
Nparams = length(x0);
% check dimensions of bounds
if( size(xbnds,1) ~= Nparams |  size(xbnds,2) ~= 2)
    disp('Dimension of xbnds is not valid')
    return
end
% TODO could check that x0 lies within bounds

x_keep=zeros(Nparams,Niter); 
L_keep=zeros(1,Niter); 

x = x0;
L = 0;

count=0;
xstep_int = xstep;
xstep = xbnds(:,2)' * xstep;
% Set step for negative values to be lower bounds:
xstep(2) = xbnds(2,1) * xstep_int; % dpHMM
xstep(6) = xbnds(6,1) * xstep_int; % dpSC

for k=1:Niter
    % generate proposal
    xprop = x + xstep.*(rand(Nparams,1)'-0.5);

    % check bounds
    if (xprop > xbnds(:,1)' & xprop < xbnds(:,2)' )

        % prediction
        dprop = fun(xprop, varargin{:});

        % likelihood
        % Lprop = prod(normpdf( data-dprop, 0, sigma));
        %Lprop = exp(-0.5*(norm(data-dprop))^2/sigma^2);
        % log likelihood
        
        % Maybe make this L1 norm instead of L2
        % Look at Kyle's distributions and see how smooth they are
        % Try MC hammer algorithm to parallelize 
        % Should be running 1M simulations

        Lprop = -0.5*(norm( (data-dprop)./sigma, 1))*1e-2;

        u=rand(1,1);

        %if (L==0 || u <= Lprop/L)
        if (L==0 || log(u) <= Lprop - L)
             count=count+1;
             x = xprop;
             L = Lprop;
        end 
    
    end
    x_keep(:,k) = x;
    L_keep(k) = L;
end

disp("Acceptance ratio: " + count/Niter)
