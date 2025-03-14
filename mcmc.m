function  [x_keep, L_keep, count, gps_l2, insar_l2] = mcmc(func,data,x0,xstep,xbnds,sigma,Niter,gps_weighting,varargin)
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
N_gps = length(varargin{3})*3;
N_insar = length(varargin{6});


x = x0;
prior_weight = 4e-1;
dprop = fun(x, varargin{:});
L_scaling = 1;%5e-2; %3e-4

L_gps = -0.5 * sum(((data(1:N_gps) - dprop(1:N_gps))./sigma(1:N_gps)).^2);
L_insar = -0.5 * sum(((data(N_gps+1:end) - dprop((N_gps+1:end)))./sigma((N_gps+1:end))).^2);
L = L_scaling*(0*L_gps*gps_weighting + L_insar);% + prior_weight*log(get_prior(x));
barLength = 25;
prevStr = '';



count=0;
xstep_int = xstep;
xstep = xbnds(:,2)' .* xstep_int;
% Set step for negative values to be lower bounds:
xstep(2) = xbnds(2,1) * xstep_int(2); % dpHMM
xstep(6) = xbnds(6,1) * xstep_int(6); % dpSC

for k=1:Niter
    % generate proposal
    xprop = x + xstep.*(rand(Nparams,1)'-0.5);
    % check bounds
    if all(xprop > xbnds(:,1)') && all(xprop < xbnds(:,2)')
        
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
        prior_prob = get_prior(xprop);


        L_gps = -0.5 * sum(((data(1:N_gps) - dprop(1:N_gps))./sigma(1:N_gps)).^2);
        L_insar = -0.5 * sum(((data(N_gps+1:end) - dprop((N_gps+1:end)))./sigma((N_gps+1:end))).^2);

        Lprop = L_scaling * (0*L_gps*gps_weighting + L_insar);% + prior_weight*log(prior_prob);

        u=rand(1,1);

        %if (L==0 || u <= Lprop/L)
        if (log(u) <= Lprop - L)
             count=count+1;
             x = xprop;
             L = Lprop;
        end 
    
    end
    x_keep(:,k) = x;
    L_keep(k) = L;

    if(mod(k, 1e3) == 0)
        % --- Update Progress Bar ---
        percent = k / Niter;
        acc_ratio = count/k;
        numEquals = round(percent * barLength);
        bar = [repmat('=', 1, numEquals) repmat(' ', 1, barLength - numEquals)];
        progressStr = sprintf('[%s] %3.1f%%', bar, acc_ratio*100);
        
        % Erase the previous progress bar
        if ~isempty(prevStr)
            fprintf(repmat('\b', 1, length(prevStr)));
        end
        
        % Print the new progress bar and store the string for the next iteration
        fprintf('%s', progressStr);
        prevStr = progressStr;
        
        drawnow;  % force the update
    end
end
fprintf('\n');
disp("Acceptance ratio: " + count/Niter)

[~, idx] = max(L_keep);
dprop = fun(x_keep(:,idx)', varargin{:});

gps_l2 = sum(((data(1:N_gps) - dprop(1:N_gps))./sigma(1:N_gps)).^2);
insar_l2 = sum(((data(N_gps+1:end) - dprop((N_gps+1:end)))./sigma((N_gps+1:end))).^2);

end

function p = get_prior(m)
    % Volume, dpHMM, vert semi-diam, horiz semi-diam, dip, dpSC
    mu = [3.9e9, -3e6, 277.01, 1621.5, 63, -10e6];
    sigma = [3e9, 1e6, 200, 1e3, 25, 0.5e7];
    p = mvnpdf(m, mu, diag(sigma.^2));
end